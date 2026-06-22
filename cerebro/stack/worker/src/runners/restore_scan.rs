//! `restore_scan` runner.
//!
//! Recovery producer for the archival-restore subsystem. `restore_drive` advances a
//! single file toward `Restored` by re-enqueueing itself on a delay; if that
//! poll chain is ever dropped — a worker crash, a Faktory restart, a job reaped from
//! the morgue — an archived file can be left stranded in `Requested`/`InProgress`
//! with nothing driving it. This scan is the safety net: it periodically finds those
//! in-flight restores and (re)starts a `restore_drive` for each.
//!
//! It is **idempotent**: `restore_drive` reads the authoritative state and the
//! server transitions are CAS-guarded, so re-driving an already-progressing file is
//! harmless (the poll chain simply continues). It runs on a slow recovery cadence
//! (hourly by default; see `seed.rs`) — prompt starts come from the restore-request
//! API path enqueuing `restore_drive` directly, while this scan only backstops drops.
//!
//! Only **archived** files are considered: while tiering is logical-only
//! nothing sets `archived = true`, so in production this scan is dormant but fully
//! exercisable via simulation. It never restarts terminal `Failed`/`Expired` files
//! (those require an explicit `restart = true`), so a permanently failing restore is
//! not re-driven in a loop.

use std::io;
use std::sync::Arc;
use std::time::Duration as StdDuration;

use async_trait::async_trait;
use faktory::{Job, JobRunner};
use serde::Deserialize;
use serde_json::json;

use cerebro_model::api::files::retention::RestoreState;

use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

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
    "lifecycle".to_string()
}

#[derive(Debug, Deserialize)]
struct RestoreScanArgs {
    /// Page size for the catalogue scan.
    #[serde(default = "default_page_limit")]
    page_limit: u32,
    /// Safety bound on pages scanned per run.
    #[serde(default = "default_max_pages")]
    max_pages: u32,
    /// Max `restore_drive` jobs to enqueue per run.
    #[serde(default = "default_budget")]
    budget: u32,
    /// Queue to enqueue the per-file `restore_drive` jobs on.
    #[serde(default = "default_queue")]
    queue: String,
}

/// Pure predicate: a file needs (re-)driving when it is archived and sitting in an
/// in-flight restore state. Terminal states (`Restored`, `Failed`, `Expired`) and
/// `NotArchived` are deliberately excluded — `Failed`/`Expired` require an explicit
/// restart, and a non-archived file has nothing to restore.
fn needs_drive(archived: bool, state: RestoreState) -> bool {
    archived && matches!(state, RestoreState::Requested | RestoreState::InProgress)
}

/// Periodic recovery producer: re-drive stranded in-flight restores.
pub struct RestoreScan {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
}

impl RestoreScan {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics) -> Self {
        Self { ctx, metrics }
    }

    async fn run_inner(&self, job: Job) -> Result<JobOutcome, WorkerError> {
        let args: RestoreScanArgs = parse_args(&job)?;
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
                if !needs_drive(file.archived, file.restore_state) {
                    continue;
                }
                if enqueued >= args.budget {
                    break 'outer;
                }
                let job_args = json!({ "file_id": file.id });
                match self
                    .ctx
                    .enqueue(
                        "restore_drive",
                        job_args,
                        &args.queue,
                        Some(StdDuration::from_secs(120)),
                    )
                    .await
                {
                    Ok(()) => enqueued += 1,
                    Err(e) => {
                        tracing::warn!(file = %file.id, "failed to enqueue restore_drive: {e}")
                    }
                }
            }

            if (count as u32) < args.page_limit {
                break; // short page = last page
            }
        }

        tracing::info!(scanned, enqueued, "restore_scan complete");
        Ok(JobOutcome::Succeeded)
    }
}

#[async_trait]
impl JobRunner for RestoreScan {
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
    fn drives_only_archived_in_flight() {
        assert!(needs_drive(true, RestoreState::Requested));
        assert!(needs_drive(true, RestoreState::InProgress));
        // terminal / not-archived -> not driven
        assert!(!needs_drive(true, RestoreState::Restored));
        assert!(!needs_drive(true, RestoreState::Failed));
        assert!(!needs_drive(true, RestoreState::Expired));
        assert!(!needs_drive(true, RestoreState::NotArchived));
        assert!(!needs_drive(false, RestoreState::Requested));
    }

    #[test]
    fn scan_args_defaults() {
        let a: RestoreScanArgs = serde_json::from_value(json!({})).unwrap();
        assert_eq!(a.page_limit, 200);
        assert_eq!(a.max_pages, 1000);
        assert_eq!(a.budget, 500);
        assert_eq!(a.queue, "lifecycle");
    }
}
