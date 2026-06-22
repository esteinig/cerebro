//! `restore_drive` runner.
//!
//! Drives the archival-restore state machine for a single file toward `Restored`,
//! polling the provider **by re-enqueueing itself with a delay** rather than parking
//! a worker on a multi-hour archival thaw (the blocking-work decision applied
//! to long waits).
//!
//! ## State machine (server-enforced)
//!
//! `transition_restore` (`POST /files/{id}/restore`) validates every edge against
//! [`RestoreState::can_transition_to`] with a compare-and-set, so concurrent drivers
//! cannot corrupt state and an idempotent re-run is a no-op. This runner only ever
//! follows valid edges from the *observed* state:
//!
//! ```text
//! NotArchived --(archived)--> Requested --> InProgress --> Restored
//!                                   \-----------------\--> Failed
//! Restored --> Expired ;  Failed/Expired --(restart)--> Requested
//! ```
//!
//! ## Provider seam
//!
//! The actual S3 Glacier `RestoreObject` / `head_object` status check is a
//! documented integration seam (a future `s3` feature, per FS-4). For dev/test,
//! `CEREBRO_RESTORE_SIMULATE_SECONDS` makes a restore "ready" that many seconds after
//! it was requested, exercising the whole machine without S3. With neither S3 nor the
//! simulation configured, the object stays pending and the poll budget eventually
//! fails it — honest about the missing integration rather than faking success.

use std::io;
use std::sync::Arc;
use std::time::Duration as StdDuration;

use async_trait::async_trait;
use chrono::{Duration, Utc};
use faktory::{Job, JobRunner};
use serde::Deserialize;
use serde_json::json;

use cerebro_client::client::CerebroClient;
use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::retention::RestoreState;
use cerebro_model::api::files::retention::StorageTier;
use cerebro_model::api::files::schema::FileRelocateSchema;

use crate::config::ArchiveSettings;
use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

fn default_poll_seconds() -> i64 {
    300
}
fn default_max_polls() -> u32 {
    288 // ~24h at a 5-minute poll interval
}
fn default_queue() -> String {
    "lifecycle".to_string()
}

#[derive(Debug, Deserialize)]
struct RestoreDriveArgs {
    /// File document id to restore.
    file_id: String,
    /// Seconds between status polls (re-enqueue delay).
    #[serde(default = "default_poll_seconds")]
    poll_seconds: i64,
    /// Maximum number of polls before the restore is marked `Failed` (timeout).
    #[serde(default = "default_max_polls")]
    max_polls: u32,
    /// Poll counter, carried across re-enqueues (internal).
    #[serde(default)]
    attempt: u32,
    /// Re-request a restore from a terminal `Failed`/`Expired` state (opt-in, so a
    /// stuck file is never restarted implicitly).
    #[serde(default)]
    restart: bool,
    /// Queue the poll chain re-enqueues onto.
    #[serde(default = "default_queue")]
    queue: String,
}

/// Provider restore status (the answer the archival provider gives us).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RestoreProgress {
    Pending,
    Ready,
    /// Reported by a real archival provider when a restore fails. The dev
    /// simulation never yields this; the `InProgress` handler maps it to the
    /// `Failed` state when S3 integration lands.
    #[allow(dead_code)]
    Failed,
}

/// Pure readiness decision (unit-testable without a `SeaweedFile`/context). With a
/// simulation window set, an object is `Ready` that many seconds after it was
/// requested; otherwise it stays `Pending` (no provider integration).
fn progress_from(
    simulate_seconds: Option<i64>,
    requested_at: Option<chrono::DateTime<Utc>>,
    now: chrono::DateTime<Utc>,
) -> RestoreProgress {
    match simulate_seconds {
        Some(secs) => match requested_at {
            Some(r) if now >= r + Duration::seconds(secs.max(0)) => RestoreProgress::Ready,
            _ => RestoreProgress::Pending,
        },
        None => RestoreProgress::Pending,
    }
}

/// Drive the restore state machine for one file.
pub struct RestoreDrive {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
    /// Cold-store settings. When `Some`, a ready restore re-materialises the
    /// object's bytes from the cold store before the file is declared retrievable.
    archive: Option<ArchiveSettings>,
}

impl RestoreDrive {
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

    /// Re-materialise an archived object from the cold store back into SeaweedFS and
    /// repoint the catalogue. Fetches + hash-verifies the cold bytes, writes
    /// a new local object, then relocates the file to a retrievable tier with
    /// `archived = false` and the new fid. Errors propagate so the caller marks the
    /// restore failed rather than declaring a non-retrievable object restored.
    async fn rematerialize(
        &self,
        api: &CerebroClient,
        file: &SeaweedFile,
        archive: &ArchiveSettings,
    ) -> anyhow::Result<()> {
        let archive_key = match &file.archive_key {
            Some(k) => k.clone(),
            None => anyhow::bail!(
                "archived file {} has no archive_key; cannot re-materialise from cold",
                file.id
            ),
        };

        // Fetch + verify + re-upload off-executor; returns the new local fid.
        let fs = self.ctx.fs()?;
        let archive_c = archive.clone();
        let name = file.name.clone();
        let effective = file.effective_identifier().to_string();
        let expected_hash = file.hash.clone();
        let outcome = self
            .ctx
            .run_blocking(move || {
                let store = archive_c
                    .open_store()
                    .map_err(|e| WorkerError::Other(format!("open cold store: {e}")))?;
                crate::archive::restore_object(
                    &fs,
                    store.as_ref(),
                    &archive_key,
                    &effective,
                    &name,
                    Some(&expected_hash),
                )
                .map_err(|e| WorkerError::Other(e.to_string()))
            })
            .await?;

        // Repoint to the new local copy and land at Warm. CAS guards on the file's
        // current (Cold) tier. `fid` is updated only when a fresh weed object was
        // written; a path-addressed file was overwritten in place. The cold backup
        // pointer (`archive_key`) is RETAINED so the restored file keeps a
        // durable cold copy that verify-repair can re-pull on a future mismatch.
        let id = file.id.clone();
        let expected_tier = file.tier;
        let api = api.clone();
        let schema = FileRelocateSchema {
            expected_tier,
            target_tier: StorageTier::Warm,
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
        Ok(())
    }

    async fn run_inner(&self, job: Job) -> Result<JobOutcome, WorkerError> {
        let args: RestoreDriveArgs = parse_args(&job)?;
        let api = self.ctx.api()?.clone();
        let id = args.file_id.clone();

        let file = {
            let api = api.clone();
            let id = id.clone();
            self.ctx
                .run_blocking(move || {
                    api.get_file(&id)
                        .map_err(|e| WorkerError::Api(e.to_string()))
                })
                .await?
        };

        match file.restore_state {
            RestoreState::Restored => {
                tracing::info!(%id, "already restored; object retrievable");
                Ok(JobOutcome::Skipped)
            }
            RestoreState::NotArchived => {
                if !file.archived {
                    tracing::info!(%id, "file is not archived; nothing to restore");
                    return Ok(JobOutcome::Skipped);
                }
                self.transition(&api, &id, RestoreState::Requested).await?;
                self.poll_again(&args, "requested archival restore").await
            }
            RestoreState::Requested => {
                // Provider has the request; advance to in-progress and start polling.
                self.transition(&api, &id, RestoreState::InProgress).await?;
                self.poll_again(&args, "restore in progress").await
            }
            RestoreState::InProgress => match self.restore_progress(&file) {
                RestoreProgress::Ready => {
                    // Re-materialise the bytes from cold before declaring the object
                    // retrievable. With no cold store configured this is a
                    // no-op (the simulation seam drives dev/test).
                    if let Some(archive) = &self.archive {
                        if let Err(e) = self.rematerialize(&api, &file, archive).await {
                            tracing::error!(%id, "restore re-materialisation failed: {e:#}");
                            self.transition(&api, &id, RestoreState::Failed).await?;
                            return Err(WorkerError::Other(format!(
                                "restore re-materialisation failed for {id}: {e}"
                            )));
                        }
                    }
                    self.transition(&api, &id, RestoreState::Restored).await?;
                    tracing::info!(%id, "restore complete; object retrievable");
                    Ok(JobOutcome::Succeeded)
                }
                RestoreProgress::Failed => {
                    self.transition(&api, &id, RestoreState::Failed).await?;
                    Err(WorkerError::Other(format!(
                        "archival restore failed at provider for {id}"
                    )))
                }
                RestoreProgress::Pending => {
                    if args.attempt + 1 >= args.max_polls {
                        tracing::warn!(%id, attempts = args.attempt + 1, "restore poll budget exhausted; marking failed");
                        self.transition(&api, &id, RestoreState::Failed).await?;
                        return Err(WorkerError::Other(format!(
                            "restore timed out after {} polls for {id}",
                            args.attempt + 1
                        )));
                    }
                    self.poll_again(&args, "restore still in progress").await
                }
            },
            RestoreState::Failed | RestoreState::Expired if args.restart => {
                self.transition(&api, &id, RestoreState::Requested).await?;
                self.poll_again(&args, "re-requesting restore").await
            }
            RestoreState::Failed => {
                tracing::warn!(%id, "restore previously failed; enqueue restore_drive with restart=true to retry");
                Ok(JobOutcome::Skipped)
            }
            RestoreState::Expired => {
                tracing::info!(%id, "restore window elapsed; enqueue restore_drive with restart=true to re-restore");
                Ok(JobOutcome::Skipped)
            }
        }
    }

    /// Ask the archival provider whether the object is restored. See the module
    /// docs: real S3 integration is a seam; `CEREBRO_RESTORE_SIMULATE_SECONDS`
    /// drives it for dev/test.
    fn restore_progress(&self, file: &SeaweedFile) -> RestoreProgress {
        let sim = self.ctx.config.restore_simulate_seconds;
        if sim.is_none() {
            tracing::warn!(
                file = %file.id,
                "no restore provider configured (set CEREBRO_RESTORE_SIMULATE_SECONDS for dev, or integrate S3); staying pending"
            );
        }
        progress_from(sim, file.restore_requested_at, Utc::now())
    }

    /// Apply a restore-state transition via the CAS-guarded server endpoint.
    async fn transition(
        &self,
        api: &CerebroClient,
        id: &str,
        target: RestoreState,
    ) -> Result<(), WorkerError> {
        let api = api.clone();
        let id = id.to_string();
        self.ctx
            .run_blocking(move || {
                api.transition_restore(&id, target)
                    .map(|_| ())
                    .map_err(|e| WorkerError::Api(e.to_string()))
            })
            .await
    }

    /// Re-enqueue this driver for the same file at `now + poll_seconds`, carrying the
    /// incremented attempt counter. `retry = 0`: the explicit poll chain is the only
    /// continuation, so Faktory's own retry never races it.
    async fn poll_again(
        &self,
        args: &RestoreDriveArgs,
        msg: &str,
    ) -> Result<JobOutcome, WorkerError> {
        let at = Utc::now() + Duration::seconds(args.poll_seconds.max(1));
        let next = json!({
            "file_id": args.file_id,
            "poll_seconds": args.poll_seconds,
            "max_polls": args.max_polls,
            "attempt": args.attempt + 1,
            "restart": args.restart,
            "queue": args.queue,
        });
        self.ctx
            .enqueue_at(
                "restore_drive",
                next,
                &args.queue,
                Some(at),
                Some(StdDuration::from_secs(120)),
                Some(0),
            )
            .await?;
        tracing::info!(file = %args.file_id, next_in_s = args.poll_seconds, attempt = args.attempt + 1, "{msg}; re-enqueued poll");
        Ok(JobOutcome::Succeeded)
    }
}

#[async_trait]
impl JobRunner for RestoreDrive {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn ago(secs: i64) -> Option<chrono::DateTime<Utc>> {
        Some(Utc::now() - Duration::seconds(secs))
    }

    #[test]
    fn simulation_ready_after_elapsed() {
        // requested 120s ago, sim window 60s -> ready
        assert_eq!(
            progress_from(Some(60), ago(120), Utc::now()),
            RestoreProgress::Ready
        );
        // requested 10s ago, sim window 60s -> still pending
        assert_eq!(
            progress_from(Some(60), ago(10), Utc::now()),
            RestoreProgress::Pending
        );
    }

    #[test]
    fn no_provider_stays_pending() {
        assert_eq!(
            progress_from(None, ago(10_000), Utc::now()),
            RestoreProgress::Pending
        );
    }

    #[test]
    fn not_yet_requested_is_pending() {
        assert_eq!(
            progress_from(Some(60), None, Utc::now()),
            RestoreProgress::Pending
        );
    }

    #[test]
    fn restore_args_defaults() {
        let a: RestoreDriveArgs = serde_json::from_value(json!({ "file_id": "F1" })).unwrap();
        assert_eq!(a.poll_seconds, 300);
        assert_eq!(a.max_polls, 288);
        assert_eq!(a.attempt, 0);
        assert!(!a.restart);
        assert_eq!(a.queue, "lifecycle");
    }
}
