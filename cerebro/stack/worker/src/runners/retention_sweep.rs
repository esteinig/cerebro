//! `retention_sweep` + `purge_reclaim` runners (S3-2b).
//!
//! The non-destructive → destructive halves of the retention lifecycle. Both are
//! **bulk** operations (the server endpoints act on every eligible file in one
//! pass, optionally scoped to a run/sample), so unlike the tier mover these runners
//! call the endpoint directly rather than fanning out per-file jobs.
//!
//! ## `retention_sweep` (quarantine)
//!
//! Calls `POST /files/expire` (`expire_files`): files past their `retain_until`
//! become `expiry_state = quarantined` with `quarantined_at = now`. **Non-destructive
//! and reversible** — nothing is deleted; the file simply enters the grace window.
//! Legally held files are never quarantined (enforced server-side). Safe to schedule
//! frequently.
//!
//! ## `purge_reclaim` (hard purge + byte reclamation)
//!
//! Calls `POST /files/purge` (`purge_files`): the server permanently removes the DB
//! records of quarantined files **past the grace window** (not legally held) and
//! returns their storage `fid`s. This runner then **reclaims the bytes** by deleting
//! each fid from Cerebro FS.
//!
//! Ordering is deliberate and safe: the server deletes the catalogue record *first*
//! and only then hands back fids, so the worst case if this worker dies mid-reclaim
//! is a **leaked byte object** (an orphan, surfaced by the failure metric for a
//! future reconcile pass) — never a live record pointing at deleted bytes. Byte
//! deletion is gated on `!dry_run`: a dry run returns the eligible fids and deletes
//! nothing (neither records nor bytes).

use std::io;
use std::sync::Arc;

use async_trait::async_trait;
use faktory::{Job, JobRunner};
use serde::Deserialize;

use cerebro_model::api::files::telemetry::{TelemetryEvent, TelemetryOp, TelemetryOutcome};

use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

// ===========================================================================
// retention_sweep
// ===========================================================================

#[derive(Debug, Deserialize)]
struct RetentionSweepArgs {
    /// Optional scope: only sweep this run.
    #[serde(default)]
    run_id: Option<String>,
    /// Optional scope: only sweep this sample.
    #[serde(default)]
    sample_id: Option<String>,
    /// Preview only — count what *would* be quarantined without changing state.
    #[serde(default)]
    dry_run: bool,
}

/// Quarantine files past retention (non-destructive).
pub struct RetentionSweep {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
}

impl RetentionSweep {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics) -> Self {
        Self { ctx, metrics }
    }

    async fn run_inner(&self, job: Job) -> Result<JobOutcome, WorkerError> {
        let RetentionSweepArgs {
            run_id,
            sample_id,
            dry_run,
        } = parse_args(&job)?;
        let api = self.ctx.api()?.clone();

        // Bulk server-side quarantine. The Expire telemetry is recorded server-side
        // (API-mediated), so it is not re-counted here.
        let (quarantined, dry) = {
            let api = api.clone();
            let run_id = run_id.clone();
            let sample_id = sample_id.clone();
            self.ctx
                .run_blocking(move || {
                    api.expire_files(run_id, sample_id, dry_run)
                        .map_err(|e| WorkerError::Api(e.to_string()))
                })
                .await?
        };

        tracing::info!(quarantined, dry_run = dry, "retention sweep complete");
        Ok(JobOutcome::Succeeded)
    }
}

#[async_trait]
impl JobRunner for RetentionSweep {
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

// ===========================================================================
// purge_reclaim
// ===========================================================================

#[derive(Debug, Deserialize)]
struct PurgeReclaimArgs {
    /// Optional scope: only purge this run.
    #[serde(default)]
    run_id: Option<String>,
    /// Optional scope: only purge this sample.
    #[serde(default)]
    sample_id: Option<String>,
    /// Preview only — return eligible fids without deleting records or bytes.
    #[serde(default)]
    dry_run: bool,
}

/// Hard-purge past-grace quarantined files and reclaim their bytes.
pub struct PurgeReclaim {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
}

impl PurgeReclaim {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics) -> Self {
        Self { ctx, metrics }
    }

    async fn run_inner(&self, job: Job) -> Result<JobOutcome, WorkerError> {
        let PurgeReclaimArgs {
            run_id,
            sample_id,
            dry_run,
        } = parse_args(&job)?;
        let api = self.ctx.api()?.clone();

        // 1. Server-side purge: deletes the catalogue records (Purge telemetry is
        //    recorded server-side) and returns the fids whose bytes we now own the
        //    job of reclaiming. In dry-run, nothing is deleted server-side.
        let (fids, dry) = {
            let api = api.clone();
            let run_id = run_id.clone();
            let sample_id = sample_id.clone();
            self.ctx
                .run_blocking(move || {
                    api.purge_files(run_id, sample_id, dry_run)
                        .map_err(|e| WorkerError::Api(e.to_string()))
                })
                .await?
        };

        if dry {
            tracing::info!(
                eligible = fids.len(),
                "purge_reclaim dry-run; no records or bytes deleted"
            );
            return Ok(JobOutcome::Skipped);
        }

        // 2. Reclaim bytes for each purged fid (best-effort, per-fid). A failure
        //    leaves an orphaned object — surfaced by the failure metric — rather
        //    than aborting the batch.
        let fs = self.ctx.fs()?.clone();
        let mut reclaimed = 0u64;
        let mut failed = 0u64;
        for fid in fids {
            let res = {
                let fs = fs.clone();
                let fid_c = fid.clone();
                self.ctx
                    .run_blocking(move || {
                        fs.delete_file(&fid_c)
                            .map_err(|e| WorkerError::Fs(e.to_string()))
                    })
                    .await
            };
            match res {
                Ok(()) => {
                    reclaimed += 1;
                    self.metrics.record(&TelemetryEvent::with_detail(
                        TelemetryOp::Delete,
                        TelemetryOutcome::Success,
                        "reclaim",
                    ));
                }
                Err(e) => {
                    failed += 1;
                    tracing::warn!(%fid, "byte reclaim failed (orphaned until reconcile): {e}");
                    self.metrics.record(&TelemetryEvent::with_detail(
                        TelemetryOp::Delete,
                        TelemetryOutcome::Failure,
                        "reclaim",
                    ));
                }
            }
        }

        tracing::info!(reclaimed, failed, "purge_reclaim complete");
        Ok(JobOutcome::Succeeded)
    }
}

#[async_trait]
impl JobRunner for PurgeReclaim {
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
    fn sweep_args_default_to_unscoped_live_run() {
        let a: RetentionSweepArgs = serde_json::from_value(serde_json::json!({})).unwrap();
        assert!(a.run_id.is_none() && a.sample_id.is_none());
        assert!(!a.dry_run, "live (non-dry) by default");
    }

    #[test]
    fn purge_args_parse_scope_and_dry_run() {
        let a: PurgeReclaimArgs =
            serde_json::from_value(serde_json::json!({ "run_id": "R1", "dry_run": true })).unwrap();
        assert_eq!(a.run_id.as_deref(), Some("R1"));
        assert!(a.sample_id.is_none());
        assert!(a.dry_run);
    }
}
