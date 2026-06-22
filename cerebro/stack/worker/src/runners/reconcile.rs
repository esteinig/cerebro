//! `reconcile_scan` + `reconcile_reclaim` (S4-3).
//!
//! * `reconcile_scan` — **report-first** (D6). Enumerates the catalogue, probes
//!   each non-archived entry's object presence, and reports dangling references
//!   (catalogue entries whose object is missing). It never mutates state; S4-5
//!   repairs from the report. Safe to run on a schedule.
//! * `reconcile_reclaim` — **operator-gated** and destructive. Deletes an explicit
//!   list of orphan objects, only when enqueued with `confirm=true`. Not seeded;
//!   an operator triggers it after reviewing a report.

use std::collections::HashSet;
use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use async_trait::async_trait;
use chrono::{DateTime, Duration, NaiveDate, Utc};
use faktory::{Job, JobRunner};
use serde::Deserialize;

use crate::backup::{FilesystemObjectStore, ObjectStore};
use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::reconcile::{
    find_dangling, find_orphans, CatalogueRef, ReconcileReport, StoreObjectRef,
};
use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

const DEFAULT_BUDGET: u32 = 5_000;
const DEFAULT_GRACE_DAYS: i64 = 1;
const PAGE_LIMIT: u32 = 500;
/// Upper bound on objects walked in one store enumeration (H2). A truncated walk
/// still yields only true orphans (each is checked against the full catalogue), so
/// this just caps the work per run.
const STORE_ENUM_BUDGET: usize = 100_000;
const DEFAULT_MAX_DELETE: usize = 100;

#[derive(Debug, Deserialize)]
struct ScanArgs {
    #[serde(default = "default_budget")]
    budget: u32,
    #[serde(default = "default_grace")]
    grace_days: i64,
}
fn default_budget() -> u32 {
    DEFAULT_BUDGET
}
fn default_grace() -> i64 {
    DEFAULT_GRACE_DAYS
}

/// Parse the catalogue `date` (e.g. `2025-01-01`) to a UTC timestamp; `None` when
/// it is not in the expected day format (then grace can't protect that entry).
fn parse_catalogue_date(s: &str) -> Option<DateTime<Utc>> {
    NaiveDate::parse_from_str(s, "%Y-%m-%d")
        .ok()
        .and_then(|d| d.and_hms_opt(0, 0, 0))
        .map(|naive| DateTime::<Utc>::from_naive_utc_and_offset(naive, Utc))
}

pub struct ReconcileScan {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
    /// Object-store root for persisting reports (reuses the S4-2 backup store when
    /// configured). `None` => reports are logged only.
    report_sink: Option<PathBuf>,
}

impl ReconcileScan {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics, report_sink: Option<PathBuf>) -> Self {
        Self {
            ctx,
            metrics,
            report_sink,
        }
    }

    async fn run_inner(&self, job: Job) -> anyhow::Result<()> {
        let args: ScanArgs = parse_args(&job).unwrap_or(ScanArgs {
            budget: DEFAULT_BUDGET,
            grace_days: DEFAULT_GRACE_DAYS,
        });
        let grace = Duration::days(args.grace_days.max(0));
        let now = Utc::now();
        let budget = args.budget;

        // Enumerate the catalogue and probe object presence (blocking clients).
        let api = self.ctx.api()?;
        let fs = self.ctx.fs()?;
        let supports_enumeration = fs.enumeration_supported();
        let fs_probe = fs.clone();
        let (catalogue_count, catalogue_complete, refs, present, catalogue_paths, probe_errors) =
            self.ctx
                .run_blocking(move || {
                    let mut refs: Vec<CatalogueRef> = Vec::new();
                    let mut present: HashSet<String> = HashSet::new();
                    let mut catalogue_paths: HashSet<String> = HashSet::new();
                    let mut probe_errors = 0usize;
                    let mut catalogue_count = 0usize;
                    let mut complete = false;
                    let mut page = 0u32;
                    loop {
                        let files = api
                            .list_files(None, None, page, PAGE_LIMIT, false)
                            .map_err(|e| WorkerError::Api(e.to_string()))?;
                        if files.is_empty() {
                            complete = true; // reached the end of the catalogue
                            break;
                        }
                        for f in &files {
                            catalogue_count += 1;
                            // Full path set (cheap) for orphan detection.
                            if let Some(p) = &f.path {
                                catalogue_paths.insert(p.clone());
                            }
                            // Dangling probing is bounded by the budget (object_exists is
                            // a network round-trip per file); path collection above is not.
                            if (catalogue_count as u32) <= budget {
                                let r = CatalogueRef {
                                    file_id: f.id.clone(),
                                    fid: f.fid.clone(),
                                    created_at: parse_catalogue_date(&f.date),
                                    archived: f.archived,
                                    tier: f.tier.to_string(),
                                };
                                // Archived objects live in remote cold, not the local
                                // store: never probe them, and let the engine skip them.
                                if f.archived {
                                    refs.push(r);
                                } else {
                                    match fs_probe.object_exists(&f.fid) {
                                        Ok(true) => {
                                            present.insert(f.fid.clone());
                                            refs.push(r);
                                        }
                                        Ok(false) => refs.push(r), // definitively absent
                                        Err(e) => {
                                            // A transient fault is never treated as missing.
                                            probe_errors += 1;
                                            tracing::debug!(
                                                "object probe error for {}: {e}",
                                                f.fid
                                            );
                                        }
                                    }
                                }
                            }
                        }
                        page += 1;
                    }
                    Ok((
                        catalogue_count,
                        complete,
                        refs,
                        present,
                        catalogue_paths,
                        probe_errors,
                    ))
                })
                .await?;

        let probed_count = refs.iter().filter(|r| !r.archived).count();
        let dangling = find_dangling(&refs, &present, grace, now);

        // Report-first: detect + record, never mutate here (S4-5 repairs).
        for d in &dangling {
            tracing::warn!(
                file_id = %d.file_id, fid = %d.fid, tier = %d.tier,
                "dangling catalogue reference — backing object missing"
            );
        }

        // Orphan detection (H2): only when the store can be enumerated (filer mode)
        // AND the catalogue path set is complete — a truncated catalogue would
        // manufacture false orphans. A truncated *store* walk is fine: it only
        // misses some orphans, never invents them (every found key is checked
        // against the full catalogue).
        let mut orphans = Vec::new();
        let mut store_enumerated = false;
        let mut store_count = None;
        if supports_enumeration {
            if catalogue_complete {
                let fs_enum = fs.clone();
                match self
                    .ctx
                    .run_blocking(move || {
                        fs_enum
                            .enumerate_objects("/", STORE_ENUM_BUDGET)
                            .map_err(|e| WorkerError::Other(e.to_string()))
                    })
                    .await
                {
                    Ok(objects) => {
                        let store_refs: Vec<StoreObjectRef> = objects
                            .into_iter()
                            .map(|o| StoreObjectRef {
                                key: o.path,
                                modified_at: o.mtime,
                            })
                            .collect();
                        store_count = Some(store_refs.len());
                        orphans = find_orphans(&store_refs, &catalogue_paths, grace, now);
                        store_enumerated = true;
                        for o in &orphans {
                            tracing::warn!(key = %o.key, "orphan object — no catalogue entry");
                        }
                    }
                    Err(e) => {
                        tracing::warn!("store enumeration failed; skipping orphan detection: {e}")
                    }
                }
            } else {
                tracing::info!(
                    "catalogue exceeds reconcile budget; skipping orphan detection (raise the budget to enable it)"
                );
            }
        }

        let report = ReconcileReport {
            generated_at: now,
            catalogue_count,
            probed_count,
            store_enumerated,
            store_count,
            probe_errors,
            dangling: dangling.clone(),
            orphans: orphans.clone(),
        };
        self.persist_report(&report)?;

        tracing::info!(
            catalogue = catalogue_count,
            probed = probed_count,
            dangling = dangling.len(),
            orphans = orphans.len(),
            store_enumerated,
            errors = probe_errors,
            "reconcile scan complete"
        );
        Ok(())
    }

    fn persist_report(&self, report: &ReconcileReport) -> anyhow::Result<()> {
        if let Some(root) = &self.report_sink {
            let store = FilesystemObjectStore::new(root);
            let key = format!(
                "reconcile/{}.json",
                report.generated_at.format("%Y%m%dT%H%M%SZ")
            );
            store.put(&key, &serde_json::to_vec_pretty(report)?)?;
            tracing::info!(%key, "reconcile report written to object store");
        }
        Ok(())
    }
}

#[async_trait]
impl JobRunner for ReconcileScan {
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
                tracing::error!("reconcile_scan failed: {e:#}");
                Err(io::Error::new(io::ErrorKind::Other, e.to_string()))
            }
        }
    }
}

#[derive(Debug, Deserialize)]
struct ReclaimArgs {
    /// Must be `true` — the destructive action is operator-gated (D6).
    #[serde(default)]
    confirm: bool,
    /// Explicit orphan object keys (fids) to delete. Supplied by the operator from
    /// a reviewed reconcile report.
    #[serde(default)]
    keys: Vec<String>,
    /// Upper bound on deletions in one run (defence in depth).
    #[serde(default = "default_max_delete")]
    max_delete: usize,
}
fn default_max_delete() -> usize {
    DEFAULT_MAX_DELETE
}

pub struct ReconcileReclaim {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
}

impl ReconcileReclaim {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics) -> Self {
        Self { ctx, metrics }
    }

    async fn run_inner(&self, job: Job) -> anyhow::Result<()> {
        let args: ReclaimArgs = parse_args(&job)?;
        if !args.confirm {
            anyhow::bail!(
                "reconcile_reclaim is destructive and operator-gated: enqueue with \
                 confirm=true and an explicit `keys` list (from a reviewed report)"
            );
        }
        if args.keys.is_empty() {
            tracing::info!("reconcile_reclaim: no keys supplied; nothing to do");
            return Ok(());
        }
        let keys: Vec<String> = args.keys.into_iter().take(args.max_delete).collect();
        let requested = keys.len();

        let fs = self.ctx.fs()?;
        let deleted = self
            .ctx
            .run_blocking(move || {
                let mut deleted = 0usize;
                for k in &keys {
                    match fs.delete_store_object(k) {
                        Ok(()) => {
                            deleted += 1;
                            tracing::warn!(fid = %k, "reclaimed orphan object (operator-confirmed)");
                        }
                        Err(e) => tracing::error!(fid = %k, "orphan reclaim failed: {e}"),
                    }
                }
                Ok::<usize, WorkerError>(deleted)
            })
            .await?;

        tracing::info!(requested, deleted, "reconcile reclaim complete");
        Ok(())
    }
}

#[async_trait]
impl JobRunner for ReconcileReclaim {
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
                tracing::error!("reconcile_reclaim failed: {e:#}");
                Err(io::Error::new(io::ErrorKind::Other, e.to_string()))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalogue_date_parses_day_format() {
        assert!(parse_catalogue_date("2025-01-01").is_some());
        assert!(parse_catalogue_date("not-a-date").is_none());
        assert!(parse_catalogue_date("").is_none());
    }
}
