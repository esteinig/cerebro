//! `verify_file` + `verify_scan` runners (S3-3a).
//!
//! Scheduled integrity verification. Mirrors the mover's producer/worker split: a
//! cheap periodic **scan** selects files and fans out per-file **verify** jobs, each
//! of which re-reads the bytes and checks them against the catalogue's BLAKE3 hash.
//!
//! ## `verify_file` (per file)
//!
//! Downloads the file by fid with `verify = true` (the same integrity primitive the
//! tier mover's optional gate uses) — a hash mismatch or a missing object is an
//! integrity failure. The job itself **succeeds** (it performed the check); the
//! *result* is recorded on the lifecycle metric
//! (`cerebro_file_lifecycle_ops_total{op="verify", outcome="success|failure"}`), so a
//! corrupt file raises a `failure` counter to alert on rather than a retry storm.
//! Archival/cold objects are skipped (verifying them would require a restore).
//!
//! ## `verify_scan` (periodic producer)
//!
//! Pages the catalogue, skips cold/archived files, and enqueues one `verify_file`
//! per directly-retrievable file up to a per-run budget. Verifying the whole estate
//! every run is intentionally avoided — schedule it to sweep a bounded batch each
//! pass.
//!
//! ## Repair
//!
//! Detection only. Automated repair (re-replicate from a healthy copy) needs real
//! replication; the interim single-node `000` topology has no second copy to heal
//! from, so a failure is surfaced (metric + error log) for operator action. Repair
//! is a Stage 4 (backup/recovery) concern.

use std::io;
use std::sync::Arc;
use std::time::Duration as StdDuration;

use async_trait::async_trait;
use faktory::{Job, JobRunner};
use serde::Deserialize;
use serde_json::json;

use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::retention::StorageTier;
use cerebro_model::api::files::telemetry::{TelemetryEvent, TelemetryOp, TelemetryOutcome};

use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

/// A cold/archived object is not directly retrievable — verifying it would force a
/// restore, so the scan skips it and `verify_file` no-ops on it.
fn directly_retrievable(file: &SeaweedFile) -> bool {
    retrievable(file.tier, file.archived)
}

/// Pure predicate (unit-testable without a `SeaweedFile`).
fn retrievable(tier: StorageTier, archived: bool) -> bool {
    tier != StorageTier::Cold && !archived
}

// ===========================================================================
// verify_file
// ===========================================================================

#[derive(Debug, Deserialize)]
struct VerifyFileArgs {
    file_id: String,
}

/// Verify a single file's bytes against its catalogue hash.
pub struct VerifyFile {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
}

impl VerifyFile {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics) -> Self {
        Self { ctx, metrics }
    }

    async fn run_inner(&self, job: Job) -> Result<JobOutcome, WorkerError> {
        let VerifyFileArgs { file_id } = parse_args(&job)?;
        let api = self.ctx.api()?.clone();

        let file = {
            let api = api.clone();
            let id = file_id.clone();
            self.ctx
                .run_blocking(move || api.get_file(&id).map_err(|e| WorkerError::Api(e.to_string())))
                .await?
        };

        if !directly_retrievable(&file) {
            tracing::info!(%file_id, "skipping verify of archived/cold file (needs restore)");
            return Ok(JobOutcome::Skipped);
        }

        // Download with verify=true: re-reads the bytes and checks BLAKE3 == stored hash.
        let fs = self.ctx.fs()?.clone();
        let fid = file.fid.clone();
        let tmp = std::env::temp_dir().join(format!("cerebro-verify-{}", file.id));
        let outdir = tmp.clone();
        let result = self
            .ctx
            .run_blocking(move || {
                std::fs::create_dir_all(&outdir)
                    .map_err(|e| WorkerError::Fs(format!("create temp dir: {e}")))?;
                fs.download_files(&vec![fid], None, None, &outdir, true)
                    .map(|_| ())
                    .map_err(|e| WorkerError::Fs(e.to_string()))
            })
            .await;
        let _ = std::fs::remove_dir_all(&tmp);

        match result {
            Ok(()) => {
                self.metrics.record(&TelemetryEvent::success(TelemetryOp::Verify));
                tracing::info!(%file_id, "integrity verified");
            }
            Err(e) => {
                // The job ran fine; the *file* failed. Surface on the lifecycle
                // metric (alert target) + error log, not as a job failure (which
                // would retry-storm on genuinely corrupt bytes).
                self.metrics.record(&TelemetryEvent::failure(TelemetryOp::Verify));
                tracing::error!(%file_id, fid = %file.fid, "INTEGRITY FAILURE on verify: {e}");
            }
        }
        Ok(JobOutcome::Succeeded)
    }
}

#[async_trait]
impl JobRunner for VerifyFile {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {
        let kind = job.kind().to_string();
        self.metrics.record_job(&kind, JobOutcome::Started);
        match self.run_inner(job).await {
            Ok(JobOutcome::Skipped) => {
                self.metrics.record_job(&kind, JobOutcome::Skipped);
                Ok(())
            }
            Ok(_) => {
                self.metrics.record_job(&kind, JobOutcome::Succeeded);
                Ok(())
            }
            Err(e) => {
                self.metrics.record_job(&kind, JobOutcome::Failed);
                Err(e.into())
            }
        }
    }
}

// ===========================================================================
// verify_scan
// ===========================================================================

fn default_page_limit() -> u32 {
    200
}
fn default_max_pages() -> u32 {
    1000
}
fn default_budget() -> u32 {
    500
}
fn default_queue() -> String {
    "maintenance".to_string()
}

#[derive(Debug, Deserialize)]
struct VerifyScanArgs {
    /// Page size for the catalogue scan.
    #[serde(default = "default_page_limit")]
    page_limit: u32,
    /// Safety bound on pages scanned per run.
    #[serde(default = "default_max_pages")]
    max_pages: u32,
    /// Max `verify_file` jobs to enqueue per run (bounds the work per sweep).
    #[serde(default = "default_budget")]
    budget: u32,
    /// Queue to enqueue the per-file verify jobs on.
    #[serde(default = "default_queue")]
    queue: String,
}

/// Periodic producer: enqueue `verify_file` for a bounded batch of directly
/// retrievable files.
pub struct VerifyScan {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
}

impl VerifyScan {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics) -> Self {
        Self { ctx, metrics }
    }

    async fn run_inner(&self, job: Job) -> Result<JobOutcome, WorkerError> {
        let args: VerifyScanArgs = parse_args(&job)?;
        let api = self.ctx.api()?.clone();

        let mut scanned = 0u64;
        let mut enqueued = 0u32;

        'outer: for page in 0..args.max_pages {
            let files = {
                let api = api.clone();
                let limit = args.page_limit;
                self.ctx
                    .run_blocking(move || {
                        api.list_files(None, None, page, limit, false)
                            .map_err(|e| WorkerError::Api(e.to_string()))
                    })
                    .await?
            };
            if files.is_empty() {
                break;
            }
            let count = files.len();
            scanned += count as u64;

            for file in &files {
                if !directly_retrievable(file) {
                    continue;
                }
                if enqueued >= args.budget {
                    break 'outer;
                }
                let job_args = json!({ "file_id": file.id });
                match self
                    .ctx
                    .enqueue("verify_file", job_args, &args.queue, Some(StdDuration::from_secs(1800)))
                    .await
                {
                    Ok(()) => enqueued += 1,
                    Err(e) => tracing::warn!(file = %file.id, "failed to enqueue verify_file: {e}"),
                }
            }

            if (count as u32) < args.page_limit {
                break; // short page = last page
            }
        }

        tracing::info!(scanned, enqueued, "verify_scan complete");
        self.metrics
            .record(&TelemetryEvent::with_detail(TelemetryOp::Verify, TelemetryOutcome::Success, "scan"));
        Ok(JobOutcome::Succeeded)
    }
}

#[async_trait]
impl JobRunner for VerifyScan {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {
        let kind = job.kind().to_string();
        self.metrics.record_job(&kind, JobOutcome::Started);
        match self.run_inner(job).await {
            Ok(_) => {
                self.metrics.record_job(&kind, JobOutcome::Succeeded);
                Ok(())
            }
            Err(e) => {
                self.metrics.record_job(&kind, JobOutcome::Failed);
                Err(e.into())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hot_warm_are_verifiable_cold_archived_are_not() {
        assert!(retrievable(StorageTier::Hot, false));
        assert!(retrievable(StorageTier::Warm, false));
        assert!(!retrievable(StorageTier::Cold, false));
        assert!(!retrievable(StorageTier::Hot, true));
    }

    #[test]
    fn scan_args_defaults() {
        let a: VerifyScanArgs = serde_json::from_value(json!({})).unwrap();
        assert_eq!(a.page_limit, 200);
        assert_eq!(a.budget, 500);
        assert_eq!(a.queue, "maintenance");
    }
}
