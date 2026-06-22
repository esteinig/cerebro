//! `archive_reclaim` (S4-6 H3): reclaim the redundant local copy of an archived
//! file once it is safe to do so (D4).
//!
//! After S4-4 archives a file, the authoritative copy lives in the cold object
//! store but a local SeaweedFS copy is kept. This runner deletes that local copy —
//! reclaiming hot/warm space — only when **all** of the following hold:
//!
//! * the file is `archived` with an `archive_key`;
//! * its archival is older than the grace window (anchored on `tier_moved_at`);
//! * the **cold copy is confirmed present** (`ObjectStore::exists` — the safety
//!   gate, so a local copy is never deleted without its replacement); and
//! * a local copy actually still exists.
//!
//! After deletion the file stays `archived`, so reconcile skips it and retrieval
//! goes through the restore path (S4-5). The stale fid is harmless and is replaced
//! on the next restore. This pass mutates only by deleting confirmed-redundant
//! local bytes; it never touches the catalogue.

use std::io;
use std::sync::Arc;

use async_trait::async_trait;
use chrono::{DateTime, Duration, Utc};
use faktory::{Job, JobRunner};
use serde::Deserialize;

use crate::config::ArchiveSettings;
use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

const DEFAULT_BUDGET: u32 = 2_000;
const DEFAULT_MAX_DELETE: usize = 200;
const PAGE_LIMIT: u32 = 500;

/// Whether a file is eligible for local-copy reclaim (H3): it is archived, has a
/// cold key, and its archival (anchored on `tier_moved_at`) is older than the grace
/// window. Pure, so it can be unit-tested without a worker or store. The
/// cold-copy-present and local-copy-present checks are done separately by the runner
/// (they require I/O); this is the cheap catalogue-only gate.
fn is_reclaim_candidate(
    archived: bool,
    archive_key: Option<&str>,
    tier_moved_at: Option<DateTime<Utc>>,
    grace: Duration,
    now: DateTime<Utc>,
) -> bool {
    archived && archive_key.is_some() && matches!(tier_moved_at, Some(t) if now - t >= grace)
}

#[derive(Debug, Deserialize)]
struct ReclaimArgs {
    /// Catalogue entries scanned per run.
    #[serde(default = "default_budget")]
    budget: u32,
    /// Local copies reclaimed per run (defence in depth).
    #[serde(default = "default_max_delete")]
    max_delete: usize,
    /// Override the configured grace (days).
    #[serde(default)]
    grace_days: Option<i64>,
}
fn default_budget() -> u32 {
    DEFAULT_BUDGET
}
fn default_max_delete() -> usize {
    DEFAULT_MAX_DELETE
}

pub struct ArchiveReclaim {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
    archive: Option<ArchiveSettings>,
}

impl ArchiveReclaim {
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

    async fn run_inner(&self, job: Job) -> anyhow::Result<()> {
        // Reclaim needs the cold store for the safety gate; without it, do nothing.
        let archive = match &self.archive {
            Some(a) => a.clone(),
            None => {
                tracing::info!("archive_reclaim: no cold store configured; nothing to do");
                return Ok(());
            }
        };
        let args: ReclaimArgs = parse_args(&job).unwrap_or(ReclaimArgs {
            budget: DEFAULT_BUDGET,
            max_delete: DEFAULT_MAX_DELETE,
            grace_days: None,
        });
        let grace = Duration::days(args.grace_days.unwrap_or(archive.local_grace_days).max(0));
        let now = Utc::now();
        let budget = args.budget;
        let max_delete = args.max_delete;

        let api = self.ctx.api()?;
        let fs = self.ctx.fs()?;

        let (scanned, reclaimed, missing_cold) = self
            .ctx
            .run_blocking(move || {
                let store = archive
                    .open_store()
                    .map_err(|e| WorkerError::Other(format!("open cold store: {e}")))?;

                let mut scanned = 0u32;
                let mut reclaimed = 0usize;
                let mut missing_cold = 0usize;
                let mut page = 0u32;

                'outer: while scanned < budget {
                    let files = api
                        .list_files(None, None, page, PAGE_LIMIT, false)
                        .map_err(|e| WorkerError::Api(e.to_string()))?;
                    if files.is_empty() {
                        break;
                    }
                    for f in &files {
                        scanned += 1;
                        if reclaimed >= max_delete {
                            break 'outer;
                        }
                        // Candidate gate (cheap, catalogue-only): archived, has a
                        // cold key, past grace. Presence checks (cold + local) follow.
                        if !is_reclaim_candidate(
                            f.archived,
                            f.archive_key.as_deref(),
                            f.tier_moved_at,
                            grace,
                            now,
                        ) {
                            continue;
                        }
                        let key = match &f.archive_key {
                            Some(k) => k,
                            None => continue, // unreachable after the gate
                        };

                        // Safety gate: the cold copy must be confirmed present.
                        match store.exists(key) {
                            Ok(true) => {}
                            Ok(false) => {
                                missing_cold += 1;
                                tracing::warn!(
                                    file = %f.id, key = %key,
                                    "archived file missing its cold copy — NOT reclaiming local copy"
                                );
                                continue;
                            }
                            Err(e) => {
                                tracing::warn!(file = %f.id, "cold existence check failed; skipping: {e}");
                                continue;
                            }
                        }

                        // Delete by filer path when present (cleans metadata + data),
                        // else by fid. Probe the same key first so an already-reclaimed
                        // file is skipped without a spurious delete error.
                        let local_key = f.path.clone().unwrap_or_else(|| f.fid.clone());
                        match fs.store_object_present(&local_key) {
                            Ok(true) => match fs.delete_store_object(&local_key) {
                                Ok(()) => {
                                    reclaimed += 1;
                                    tracing::info!(
                                        file = %f.id, key = %local_key,
                                        "reclaimed local copy of archived file"
                                    );
                                }
                                Err(e) => tracing::warn!(file = %f.id, "local reclaim failed: {e}"),
                            },
                            Ok(false) => {} // already reclaimed
                            Err(e) => tracing::warn!(file = %f.id, "local presence check failed; skipping: {e}"),
                        }
                    }
                    page += 1;
                }
                Ok::<(u32, usize, usize), WorkerError>((scanned, reclaimed, missing_cold))
            })
            .await?;

        tracing::info!(
            scanned,
            reclaimed,
            missing_cold,
            grace_days = grace.num_days(),
            "archive reclaim complete"
        );
        Ok(())
    }
}

#[async_trait]
impl JobRunner for ArchiveReclaim {
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
                tracing::error!("archive_reclaim failed: {e:#}");
                Err(io::Error::new(io::ErrorKind::Other, e.to_string()))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_requires_archived_key_and_past_grace() {
        let now = Utc::now();
        let grace = Duration::days(7);
        // archived + key + archived 10 days ago (> 7) -> candidate
        assert!(is_reclaim_candidate(
            true,
            Some("k"),
            Some(now - Duration::days(10)),
            grace,
            now
        ));
        // still within grace (3 days) -> not yet
        assert!(!is_reclaim_candidate(
            true,
            Some("k"),
            Some(now - Duration::days(3)),
            grace,
            now
        ));
        // not archived -> never
        assert!(!is_reclaim_candidate(
            false,
            Some("k"),
            Some(now - Duration::days(10)),
            grace,
            now
        ));
        // no archive_key -> nothing to fall back to, never
        assert!(!is_reclaim_candidate(
            true,
            None,
            Some(now - Duration::days(10)),
            grace,
            now
        ));
        // unknown archival time -> cannot age it out, never
        assert!(!is_reclaim_candidate(true, Some("k"), None, grace, now));
    }

    #[test]
    fn grace_boundary_is_inclusive() {
        let now = Utc::now();
        let grace = Duration::days(7);
        // archived exactly `grace` ago: now - t == grace, and the gate is `>= grace`.
        assert!(is_reclaim_candidate(
            true,
            Some("k"),
            Some(now - Duration::days(7)),
            grace,
            now
        ));
    }

    #[test]
    fn zero_grace_reclaims_immediately() {
        let now = Utc::now();
        let grace = Duration::days(0);
        assert!(is_reclaim_candidate(true, Some("k"), Some(now), grace, now));
    }
}
