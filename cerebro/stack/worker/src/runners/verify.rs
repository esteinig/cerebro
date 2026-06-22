//! `verify_file` + `verify_scan` runners (S3-3a, hardened in S3-5 #3/#6).
//!
//! Scheduled integrity verification. Mirrors the mover's producer/worker split: a
//! cheap periodic **scan** selects files and fans out per-file **verify** jobs, each
//! of which re-reads the bytes and checks them against the catalogue's BLAKE3 hash.
//!
//! ## `verify_file` (per file)
//!
//! Streams the object's bytes through BLAKE3 and compares to the stored hash
//! (`FileSystemClient::verify_object`, S3-5 #6) — **no temp-disk copy**, constant
//! memory, suitable for multi-gigabyte read sets. Outcomes:
//! * **match** — stamp `verified_at = now` (so the scan can rotate coverage) and
//!   record `op="verify", outcome="success"`.
//! * **mismatch** — a real integrity failure: record `outcome="failure"` (an alert
//!   target) and **do not** stamp `verified_at` (the file stays at the head of the
//!   rotation). The *job* still succeeds — a corrupt object must not retry-storm.
//! * **transport/IO error** — the check didn't complete; return `Err` so Faktory
//!   retries.
//! Archived objects are skipped (verifying them would require a restore). While
//! tiering is logical-only (S3-5 #1) every non-archived file is directly
//! retrievable regardless of tier, so cold logical files are still verified.
//!
//! ## `verify_scan` (periodic producer)
//!
//! Pages the catalogue and enqueues one `verify_file` per non-archived file that
//! is **due** — never verified, or last verified longer ago than
//! `min_reverify_secs` — up to a per-run budget. Because a successful verify stamps
//! `verified_at`, recently-checked files drop out of the due set and each run picks
//! up the *least-recently-verified* files, sweeping the whole estate over time
//! instead of re-checking the head of the list every run.
//!
//! ## Repair (S4-6 H4)
//!
//! When a verify mismatch hits a file that has a cold backup (`archive_key`, now
//! retained across restores), `verify_file` repairs it: it re-pulls the bytes from
//! cold, verifies them against the catalogue hash, writes them back to the file's
//! effective location (overwriting the filer object in place, or writing a fresh
//! weed fid), and repoints the catalogue — keeping the file live and its cold
//! backup retained. A file with no cold backup, or a repair that errors, is
//! surfaced (metric + error log) for operator action as before.

use std::io;
use std::sync::Arc;
use std::time::Duration as StdDuration;

use async_trait::async_trait;
use chrono::{DateTime, Utc};
use faktory::{Job, JobRunner};
use serde::Deserialize;
use serde_json::json;

use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::retention::StorageTier;
use cerebro_model::api::files::schema::FileRelocateSchema;
use cerebro_model::api::files::telemetry::{TelemetryEvent, TelemetryOp, TelemetryOutcome};

use crate::config::ArchiveSettings;
use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

/// Default re-verify interval: a file verified within this window is not re-checked
/// this pass (7 days). Drives rotation together with the per-run budget.
fn default_min_reverify_secs() -> i64 {
    7 * 24 * 3600
}

/// An archived object is not directly retrievable — verifying it would force a
/// restore, so the scan skips it and `verify_file` no-ops on it. With logical-only
/// tiering (S3-5 #1) the storage tier does not affect retrievability; only the real
/// `archived` flag does.
fn directly_retrievable(file: &SeaweedFile) -> bool {
    retrievable(file.archived)
}

/// Pure predicate (unit-testable without a `SeaweedFile`).
fn retrievable(archived: bool) -> bool {
    !archived
}

/// Whether a verify mismatch can be repaired from cold (H4): a cold store is
/// configured AND the file has a retained cold backup (`archive_key`). Pure, so the
/// gate is unit-testable without a worker or store.
fn is_repairable(archive_configured: bool, archive_key: Option<&str>) -> bool {
    archive_configured && archive_key.is_some()
}

/// Whether a file is due for (re-)verification: never verified, or last verified
/// longer ago than `min_interval_secs`. Pure for unit testing.
fn is_verify_due(
    verified_at: Option<DateTime<Utc>>,
    now: DateTime<Utc>,
    min_interval_secs: i64,
) -> bool {
    match verified_at {
        None => true,
        Some(ts) => (now - ts).num_seconds() >= min_interval_secs,
    }
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
    /// Cold store config (S4-6 H4). When set, a verify mismatch on a file that has
    /// a cold backup (`archive_key`) is repaired from cold instead of only alerting.
    archive: Option<ArchiveSettings>,
}

impl VerifyFile {
    pub fn new(
        ctx: Arc<WorkerContext>,
        metrics: Metrics,
        archive: Option<ArchiveSettings>,
    ) -> Self {
        Self {
            ctx,
            metrics,
            archive,
        }
    }

    /// Attempt to repair a corrupt file from its cold backup (S4-6 H4).
    ///
    /// Re-fetches the bytes from cold, verifies them against the catalogue hash
    /// (so a corrupt cold copy can never overwrite the live entry), writes them
    /// back to the file's effective location, and repoints the catalogue — keeping
    /// the file live (`archived = false`), on its current tier, with its cold backup
    /// retained. Returns `Ok(true)` on a successful repair, `Ok(false)` if no repair
    /// was possible (no cold backup configured / no `archive_key`), or `Err` if the
    /// repair was attempted but failed.
    async fn repair_from_cold(&self, file: &SeaweedFile) -> anyhow::Result<bool> {
        if !is_repairable(self.archive.is_some(), file.archive_key.as_deref()) {
            return Ok(false); // no cold store configured, or no retained backup
        }
        let archive = self
            .archive
            .clone()
            .expect("repairable checked archive present");
        let archive_key = file
            .archive_key
            .clone()
            .expect("repairable checked key present");

        let fs = self.ctx.fs()?;
        let name = file.name.clone();
        let effective = file.effective_identifier().to_string();
        let expected_hash = file.hash.clone();
        let key = archive_key.clone();
        let outcome = self
            .ctx
            .run_blocking(move || {
                let store = archive
                    .open_store()
                    .map_err(|e| WorkerError::Other(format!("open cold store: {e}")))?;
                crate::archive::restore_object(
                    &fs,
                    store.as_ref(),
                    &key,
                    &effective,
                    &name,
                    Some(&expected_hash),
                )
                .map_err(|e| WorkerError::Other(e.to_string()))
            })
            .await?;

        // Repoint to the repaired copy. Stay live, stay on the current tier, RETAIN
        // the cold backup. `fid` is updated only when a fresh weed object was written
        // (a path-addressed file was overwritten in place). CAS guards on the tier.
        let api = self.ctx.api()?;
        let id = file.id.clone();
        let wrote_new_fid = outcome.new_fid.is_some();
        let schema = FileRelocateSchema {
            expected_tier: file.tier,
            target_tier: file.tier,
            archived: false,
            fid: outcome.new_fid,
            archive_key: None,
            clear_archive_key: false,
        };
        self.ctx
            .run_blocking(move || {
                api.relocate_file(&id, &schema)
                    .map_err(|e| WorkerError::Api(e.to_string()))
            })
            .await?;

        // A fresh weed object replaced the corrupt one — delete the now-orphaned old
        // object (known-bad bytes). Best-effort; an in-place filer overwrite leaves no
        // orphan, so this only runs in the weed case. Done after the repoint commits.
        if wrote_new_fid {
            let fs2 = self.ctx.fs()?;
            let old = file.fid.clone();
            if let Err(e) = self
                .ctx
                .run_blocking(move || {
                    fs2.delete_store_object(&old)
                        .map_err(|e| WorkerError::Fs(e.to_string()))
                })
                .await
            {
                tracing::warn!(file = %file.id, "repaired but failed to delete old corrupt object: {e}");
            }
        }
        Ok(true)
    }

    async fn run_inner(&self, job: Job) -> Result<JobOutcome, WorkerError> {
        let VerifyFileArgs { file_id } = parse_args(&job)?;
        let api = self.ctx.api()?.clone();

        let file = {
            let api = api.clone();
            let id = file_id.clone();
            self.ctx
                .run_blocking(move || {
                    api.get_file(&id)
                        .map_err(|e| WorkerError::Api(e.to_string()))
                })
                .await?
        };

        if !directly_retrievable(&file) {
            tracing::info!(%file_id, "skipping verify of archived file (needs restore)");
            return Ok(JobOutcome::Skipped);
        }

        // Stream the stored bytes through BLAKE3 (no temp disk) and compare.
        let fs = self.ctx.fs()?.clone();
        let identifier = file.effective_identifier().to_string();
        let expected = file.hash.clone();
        let verify_result = self
            .ctx
            .run_blocking(move || {
                fs.verify_object(&identifier, &expected)
                    .map_err(|e| WorkerError::Fs(e.to_string()))
            })
            .await;

        match verify_result {
            Ok(true) => {
                // Stamp verified_at so the scan rotates to other files next pass.
                let stamp = {
                    let api = api.clone();
                    let id = file_id.clone();
                    self.ctx
                        .run_blocking(move || {
                            api.mark_verified(&id)
                                .map_err(|e| WorkerError::Api(e.to_string()))
                        })
                        .await
                };
                if let Err(e) = stamp {
                    // The bytes are good; only the bookkeeping stamp failed. Log and
                    // still count success — the file will simply be re-picked sooner.
                    tracing::warn!(%file_id, "verify ok but failed to stamp verified_at: {e}");
                }
                self.metrics
                    .record(&TelemetryEvent::success(TelemetryOp::Verify));
                tracing::info!(%file_id, "integrity verified");
                Ok(JobOutcome::Succeeded)
            }
            Ok(false) => {
                // Integrity failure. Before alerting, try to repair from cold (H4):
                // if the file has a cold backup, re-pull the verified-good bytes and
                // repoint. Only a real failure (no backup, or repair errored) is
                // surfaced on the lifecycle metric (alert target).
                match self.repair_from_cold(&file).await {
                    Ok(true) => {
                        let stamp = {
                            let api = api.clone();
                            let id = file_id.clone();
                            self.ctx
                                .run_blocking(move || {
                                    api.mark_verified(&id)
                                        .map_err(|e| WorkerError::Api(e.to_string()))
                                })
                                .await
                        };
                        if let Err(e) = stamp {
                            tracing::warn!(%file_id, "repaired but failed to stamp verified_at: {e}");
                        }
                        self.metrics
                            .record(&TelemetryEvent::success(TelemetryOp::Verify));
                        tracing::warn!(%file_id, fid = %file.fid, "integrity failure REPAIRED from cold backup");
                        Ok(JobOutcome::Succeeded)
                    }
                    Ok(false) => {
                        self.metrics
                            .record(&TelemetryEvent::failure(TelemetryOp::Verify));
                        tracing::error!(%file_id, fid = %file.fid, "INTEGRITY FAILURE on verify: stored bytes do not match catalogue hash (no cold backup to repair from)");
                        Ok(JobOutcome::Succeeded)
                    }
                    Err(e) => {
                        self.metrics
                            .record(&TelemetryEvent::failure(TelemetryOp::Verify));
                        tracing::error!(%file_id, fid = %file.fid, "INTEGRITY FAILURE on verify: repair from cold failed: {e:#}");
                        Ok(JobOutcome::Succeeded)
                    }
                }
            }
            Err(e) => {
                // Transport/IO error — the check did not complete. Retry.
                Err(e)
            }
        }
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
    /// Re-verify interval: skip files verified more recently than this (seconds).
    #[serde(default = "default_min_reverify_secs")]
    min_reverify_secs: i64,
    /// Queue to enqueue the per-file verify jobs on.
    #[serde(default = "default_queue")]
    queue: String,
}

/// Periodic producer: enqueue `verify_file` for a bounded batch of the
/// least-recently-verified, directly-retrievable files.
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
        let now = Utc::now();

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
                // Rotation: only (re-)verify files that are due. Recently-verified
                // files are skipped, so each run sweeps different (older) files.
                if !is_verify_due(file.verified_at, now, args.min_reverify_secs) {
                    continue;
                }
                if enqueued >= args.budget {
                    break 'outer;
                }
                let job_args = json!({ "file_id": file.id });
                match self
                    .ctx
                    .enqueue(
                        "verify_file",
                        job_args,
                        &args.queue,
                        Some(StdDuration::from_secs(1800)),
                    )
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
        self.metrics.record(&TelemetryEvent::with_detail(
            TelemetryOp::Verify,
            TelemetryOutcome::Success,
            "scan",
        ));
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
    fn only_archived_is_unverifiable_under_logical_tiering() {
        assert!(retrievable(false));
        assert!(!retrievable(true));
    }

    #[test]
    fn repairable_requires_cold_store_and_retained_backup() {
        assert!(is_repairable(true, Some("archive/team/f"))); // both present
        assert!(!is_repairable(false, Some("archive/team/f"))); // no cold store
        assert!(!is_repairable(true, None)); // no retained backup
        assert!(!is_repairable(false, None));
    }

    #[test]
    fn verify_due_when_never_or_stale() {
        let now = Utc::now();
        // Never verified -> due.
        assert!(is_verify_due(None, now, 3600));
        // Verified just now -> not due.
        assert!(!is_verify_due(Some(now), now, 3600));
        // Verified 2h ago, interval 1h -> due.
        assert!(is_verify_due(
            Some(now - chrono::Duration::hours(2)),
            now,
            3600
        ));
        // Verified 30m ago, interval 1h -> not due.
        assert!(!is_verify_due(
            Some(now - chrono::Duration::minutes(30)),
            now,
            3600
        ));
    }

    #[test]
    fn scan_args_defaults() {
        let a: VerifyScanArgs = serde_json::from_value(json!({})).unwrap();
        assert_eq!(a.page_limit, 200);
        assert_eq!(a.budget, 500);
        assert_eq!(a.min_reverify_secs, 7 * 24 * 3600);
        assert_eq!(a.queue, "maintenance");
    }
}
