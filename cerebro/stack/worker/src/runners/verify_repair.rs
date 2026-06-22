//! `verify_repair`: consume the latest reconcile report and act on its
//! dangling references.
//!
//! Reconcile detects, report-first; this pass **confirms and escalates**. For each
//! dangling reference it re-checks the live state, which collapses to three cases:
//!
//! * **Recovered** — the object reappeared since the scan (a transient miss, or
//!   SeaweedFS volume replication healed it). No action.
//! * **Archived since** — the file was tier-moved to cold after the scan, so its
//!   object is now *correctly* absent locally. No action (the cold copy is
//!   authoritative; on-demand restore is the `restore_drive` loop).
//! * **Confirmed loss** — still non-archived and still missing. There is no
//!   in-system recovery source (the SeaweedFS replica shares the fid, and there is
//!   no object-level backup), so this is escalated loudly for the operator.
//!
//! Restoring *archived* objects from cold (the active repair) is the `restore_drive`
//! loop; restoring lost *redundancy* is handled by SeaweedFS volume replication.
//! This runner is the escalation pass that turns a passive weekly report into a
//! confirmed, transient-filtered data-loss signal. It mutates nothing.

use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use async_trait::async_trait;
use faktory::{Job, JobRunner};

use cerebro_client::client::CerebroClient;
use cerebro_fs::client::FileSystemClient;
use cerebro_model::api::files::model::SeaweedFile;

use crate::backup::{FilesystemObjectStore, ObjectStore};
use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::reconcile::ReconcileReport;
use crate::telemetry::{JobOutcome, Metrics};

pub struct VerifyRepair {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
    /// Object-store root where reconcile reports are persisted (reuses the backup
    /// store, as `reconcile_scan` does). `None` => nothing to consume.
    report_sink: Option<PathBuf>,
}

impl VerifyRepair {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics, report_sink: Option<PathBuf>) -> Self {
        Self {
            ctx,
            metrics,
            report_sink,
        }
    }

    async fn run_inner(&self, _job: Job) -> anyhow::Result<()> {
        let root = match &self.report_sink {
            Some(r) => r.clone(),
            None => {
                tracing::info!("verify_repair: no report store configured; nothing to consume");
                return Ok(());
            }
        };

        // Load the most recent reconcile report (keys are timestamped, so the
        // lexicographically-greatest is newest).
        let loaded = self
            .ctx
            .run_blocking(move || {
                let store = FilesystemObjectStore::new(root);
                let mut keys = store
                    .list("reconcile/")
                    .map_err(|e| WorkerError::Other(e.to_string()))?;
                keys.sort();
                let latest = match keys.last() {
                    Some(k) => k.clone(),
                    None => return Ok::<Option<(String, ReconcileReport)>, WorkerError>(None),
                };
                let bytes = store
                    .get(&latest)
                    .map_err(|e| WorkerError::Other(e.to_string()))?;
                let report: ReconcileReport = serde_json::from_slice(&bytes).map_err(|e| {
                    WorkerError::Other(format!("parse reconcile report {latest}: {e}"))
                })?;
                Ok(Some((latest, report)))
            })
            .await?;

        let (report_key, report) = match loaded {
            Some(r) => r,
            None => {
                tracing::info!("verify_repair: no reconcile reports found");
                return Ok(());
            }
        };

        if report.dangling.is_empty() {
            tracing::info!(%report_key, "verify_repair: report has no dangling references");
            return Ok(());
        }

        let api = self.ctx.api()?;
        let fs = self.ctx.fs()?;

        let mut recovered = 0u32;
        let mut archived_since = 0u32;
        let mut confirmed_loss = 0u32;

        for d in &report.dangling {
            let file = match self.get_file(&api, &d.file_id).await {
                Ok(f) => f,
                Err(e) => {
                    // Deleted or transiently unreachable since the scan — don't escalate.
                    tracing::warn!(file = %d.file_id, "verify_repair: get_file failed; skipping: {e}");
                    continue;
                }
            };

            if file.archived {
                archived_since += 1;
                tracing::info!(file = %file.id, "verify_repair: archived since scan; object correctly in cold");
                continue;
            }

            match self.object_present(&fs, &file.fid).await {
                Ok(true) => {
                    recovered += 1;
                    tracing::info!(file = %file.id, "verify_repair: object reappeared since scan; no action");
                }
                Ok(false) => {
                    confirmed_loss += 1;
                    tracing::error!(
                        file = %file.id, fid = %file.fid, tier = %d.tier,
                        "verify_repair: CONFIRMED data loss — object missing with no in-system \
                         recovery source; operator action / external object backup required"
                    );
                }
                Err(e) => {
                    // Inconclusive probe: never escalate on a transient fault.
                    tracing::warn!(file = %file.id, "verify_repair: re-probe failed; not escalating: {e}");
                }
            }
        }

        tracing::info!(
            %report_key,
            dangling = report.dangling.len(),
            recovered,
            archived_since,
            confirmed_loss,
            "verify_repair complete"
        );
        Ok(())
    }

    async fn get_file(&self, api: &CerebroClient, id: &str) -> Result<SeaweedFile, WorkerError> {
        let api = api.clone();
        let id = id.to_string();
        self.ctx
            .run_blocking(move || {
                api.get_file(&id)
                    .map_err(|e| WorkerError::Api(e.to_string()))
            })
            .await
    }

    async fn object_present(&self, fs: &FileSystemClient, fid: &str) -> Result<bool, WorkerError> {
        let fs = fs.clone();
        let fid = fid.to_string();
        self.ctx
            .run_blocking(move || {
                fs.object_exists(&fid)
                    .map_err(|e| WorkerError::Other(e.to_string()))
            })
            .await
    }
}

#[async_trait]
impl JobRunner for VerifyRepair {
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
                tracing::error!("verify_repair failed: {e:#}");
                Err(io::Error::new(io::ErrorKind::Other, e.to_string()))
            }
        }
    }
}
