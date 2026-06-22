//! `tier_move` + `tier_move_scan` runners.
//!
//! # Tier is a logical placement label
//!
//! `tier` is a **logical** placement label on the file document. The lifecycle API
//! flips it fid-preserving — [`UpdateFileLifecycleSchema`] carries no `fid`/`path`,
//! so there is no per-file byte repoint. The mover's job is therefore the **safe,
//! integrity-gated, two-phase transition** of that label, not byte shuffling.
//!
//! While tiering is **logical-only**, changing `tier` records intent and
//! drives retention/restore bookkeeping but does not relocate bytes — every
//! non-archived object stays directly retrievable regardless of tier. Real physical
//! tiering (moving cold data to S3/Glacier, with the `archived` flag and a fid
//! repoint) is what flips `archived = true` to make
//! the restore subsystem live. See the package README for the rationale.
//!
//! ## `tier_move` (per file)
//!
//! claim (`pending_tier` CAS on the current tier) → optional integrity gate
//! (download + BLAKE3 verify; heavy, off by default) → commit (`tier` CAS, which
//! clears the claim and stamps `tier_moved_at`) → rollback (clear the claim) on any
//! post-claim failure. **Idempotent**: a file already at the target with no claim
//! is a no-op (`Skipped`).
//!
//! ## `tier_move_scan` (periodic producer)
//!
//! Page the file catalogue, find files **due** for a warm→cold move (`tier == Warm`
//! and `reported_at + warm_days <= now`, not legally held, no claim in flight) —
//! recomputing the due time from `reported_at` since `cold_move_at` is not persisted
//! — and fan out one `tier_move` job per file. Producer only: it never moves bytes
//! itself, so a scan is cheap and safe to run on a schedule.

use std::io;
use std::sync::Arc;
use std::time::Duration as StdDuration;

use async_trait::async_trait;
use chrono::{DateTime, Duration, Utc};
use faktory::{Job, JobRunner};
use serde::Deserialize;
use serde_json::json;

use cerebro_client::client::CerebroClient;
use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::retention::StorageTier;
use cerebro_model::api::files::schema::FileRelocateSchema;
use cerebro_model::api::files::schema::UpdateFileLifecycleSchema;
use cerebro_model::api::files::telemetry::{TelemetryEvent, TelemetryOp, TelemetryOutcome};

use crate::archive::archive_key_for;
use crate::config::ArchiveSettings;
use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

/// Low-cardinality label for the destination tier (a metric `detail` value).
fn tier_label(t: StorageTier) -> &'static str {
    match t {
        StorageTier::Hot => "hot",
        StorageTier::Warm => "warm",
        StorageTier::Cold => "cold",
    }
}

// ===========================================================================
// tier_move
// ===========================================================================

#[derive(Debug, Deserialize)]
struct TierMoveArgs {
    /// File document id to move.
    file_id: String,
    /// Destination tier.
    target: StorageTier,
    /// Force the deep integrity gate for this job even when the worker default
    /// (`CEREBRO_WORKER_VERIFY_ON_MOVE`) is off.
    #[serde(default)]
    verify: bool,
}

/// Move a single file to a target tier via the safe two-phase transition.
pub struct TierMove {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
    /// Cold-store archival settings. When `Some`, a successful move to the
    /// Cold tier physically archives the object and repoints the catalogue.
    archive: Option<ArchiveSettings>,
}

impl TierMove {
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

    /// Physically archive a Cold-committed object to the cold store and repoint the
    /// catalogue. Non-fatal on failure: the file stays Cold-committed and
    /// directly retrievable, and the next move re-attempts the archival.
    async fn archive_to_cold(
        &self,
        api: &CerebroClient,
        file: &SeaweedFile,
        archive: &ArchiveSettings,
    ) -> anyhow::Result<()> {
        if file.archived {
            return Ok(()); // already archived (idempotent)
        }
        let team = self.ctx.config.team.clone().unwrap_or_default();
        let key = archive_key_for(&archive.prefix, &team, &file.id);

        // Copy bytes off-executor: open the cold store, read from SeaweedFS, put.
        let fs = self.ctx.fs()?;
        let archive_c = archive.clone();
        let fid = file.fid.clone();
        let key_c = key.clone();
        self.ctx
            .run_blocking(move || {
                let store = archive_c
                    .open_store()
                    .map_err(|e| WorkerError::Other(format!("open cold store: {e}")))?;
                crate::archive::archive_object(&fs, store.as_ref(), &fid, &key_c)
                    .map_err(|e| WorkerError::Other(e.to_string()))?;
                Ok::<(), WorkerError>(())
            })
            .await?;

        // Repoint the catalogue via the dedicated relocate endpoint: mark
        // archived and record the key. Tier stays Cold (CAS guards on it).
        let id = file.id.clone();
        let api = api.clone();
        let schema = FileRelocateSchema {
            expected_tier: StorageTier::Cold,
            target_tier: StorageTier::Cold,
            archived: true,
            fid: None,
            archive_key: Some(key),
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
        let args: TierMoveArgs = parse_args(&job)?;
        let api = self.ctx.api()?.clone();
        let id = args.file_id.clone();
        let target = args.target;

        // 1. Current state: provides both the idempotency check and the CAS
        //    `expected` tier so a concurrent mover can never double-apply.
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
        let current = file.tier;

        if current == target && file.pending_tier.is_none() {
            tracing::info!(%id, ?current, "file already at target tier; no-op");
            return Ok(JobOutcome::Skipped);
        }

        // 2. Claim: set pending_tier with a compare-and-set on the current tier.
        //    A rejected CAS means another mover won the race or the state moved —
        //    benign, recorded as `rejected` (not a failure).
        {
            let api_c = api.clone();
            let id_c = id.clone();
            self.ctx
                .run_blocking(move || {
                    api_c
                        .claim_tier_move(&id_c, target, current)
                        .map_err(|e| WorkerError::Api(e.to_string()))
                })
                .await
                .map_err(|e| {
                    tracing::warn!(%id, "tier-move claim rejected/failed: {e}");
                    self.metrics
                        .record(&TelemetryEvent::rejected(TelemetryOp::TierMove));
                    e
                })?;
        }

        // 3. Integrity gate (optional, heavy): prove the bytes are retrievable and
        //    intact before committing the new label. Off by default — deep
        //    verification is the scheduled verify worker's job.
        if args.verify || self.ctx.config.verify_on_move {
            if let Err(e) = self.verify_retrievable(&file).await {
                tracing::warn!(%id, "integrity gate failed; rolling back claim: {e}");
                self.rollback(&api, &id, current).await;
                self.metrics.record(&TelemetryEvent::with_detail(
                    TelemetryOp::TierMove,
                    TelemetryOutcome::Failure,
                    tier_label(target),
                ));
                return Err(e);
            }
        }

        // 4. Commit: set tier with the same CAS; clears the claim and stamps
        //    tier_moved_at. On any failure, roll the claim back.
        let commit = {
            let api = api.clone();
            let id = id.clone();
            self.ctx
                .run_blocking(move || {
                    api.commit_tier_move(&id, target, current)
                        .map_err(|e| WorkerError::Api(e.to_string()))
                })
                .await
        };
        if let Err(e) = commit {
            tracing::warn!(%id, "tier-move commit failed; rolling back claim: {e}");
            self.rollback(&api, &id, current).await;
            self.metrics.record(&TelemetryEvent::with_detail(
                TelemetryOp::TierMove,
                TelemetryOutcome::Failure,
                tier_label(target),
            ));
            return Err(e);
        }

        tracing::info!(%id, from = ?current, to = ?target, "tier move committed");
        self.metrics.record(&TelemetryEvent::with_detail(
            TelemetryOp::TierMove,
            TelemetryOutcome::Success,
            tier_label(target),
        ));

        // 5. Physical archival: a move to Cold with a cold store configured
        //    copies the bytes to the object store and repoints the catalogue
        //    (archived = true). Non-fatal — the tier is already committed, so an
        //    archival failure just leaves the file Cold-but-not-archived for the
        //    next move to retry.
        if target == StorageTier::Cold {
            if let Some(archive) = &self.archive {
                if let Err(e) = self.archive_to_cold(&api, &file, archive).await {
                    tracing::warn!(%id, "tier committed to Cold but archival failed (will retry): {e:#}");
                }
            }
        }

        Ok(JobOutcome::Succeeded)
    }

    /// Best-effort: clear an in-flight claim after a post-claim failure so a
    /// crashed/aborted move never strands a `pending_tier`. The CAS on
    /// `expected_tier` ensures we only clear *our* claim.
    async fn rollback(&self, api: &CerebroClient, id: &str, current: StorageTier) {
        let schema = UpdateFileLifecycleSchema {
            clear_pending_tier: true,
            clear_pending_since: true,
            expected_tier: Some(current),
            ..Default::default()
        };
        let api = api.clone();
        let id_owned = id.to_string();
        if let Err(e) = self
            .ctx
            .run_blocking(move || {
                api.update_file_lifecycle(&id_owned, &schema)
                    .map_err(|e| WorkerError::Api(e.to_string()))
            })
            .await
        {
            tracing::warn!(%id, "rollback (clear pending_tier) failed: {e}");
        }
    }

    /// Deep integrity gate: stream the file's stored bytes through BLAKE3 and
    /// compare to the catalogue hash (`FileSystemClient::verify_object`) —
    /// no temp-disk copy. Heavy for large artefacts only in transfer, not disk;
    /// off by default (the scheduled verify worker owns routine checks).
    async fn verify_retrievable(&self, file: &SeaweedFile) -> Result<(), WorkerError> {
        let fs = self.ctx.fs()?.clone();
        let identifier = file.effective_identifier().to_string();
        let expected = file.hash.clone();
        let ok = self
            .ctx
            .run_blocking(move || {
                fs.verify_object(&identifier, &expected)
                    .map_err(|e| WorkerError::Fs(e.to_string()))
            })
            .await?;
        if ok {
            Ok(())
        } else {
            Err(WorkerError::Fs(format!(
                "integrity gate failed: stored bytes for {} do not match catalogue hash",
                file.id
            )))
        }
    }
}

#[async_trait]
impl JobRunner for TierMove {
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
// tier_move_scan
// ===========================================================================

fn default_target() -> StorageTier {
    StorageTier::Cold
}
fn default_warm_days() -> i64 {
    365
}
fn default_page_limit() -> u32 {
    200
}
fn default_max_pages() -> u32 {
    1000
}
/// Default age after which an in-flight `pending_tier` claim is considered stale
/// and re-driven (1 hour). A crashed mover (died after claim, before commit or
/// rollback) would otherwise strand `pending_tier` forever.
fn default_stale_claim_secs() -> i64 {
    3600
}

#[derive(Debug, Deserialize)]
struct TierMoveScanArgs {
    /// Tier to move due files to (default `Cold`).
    #[serde(default = "default_target")]
    target: StorageTier,
    /// Warm dwell before the cold move, in days — recomputes the (unpersisted)
    /// `cold_move_at` from `reported_at`. Should match the deployment's
    /// `RetentionPolicy.warm_days`.
    #[serde(default = "default_warm_days")]
    warm_days: i64,
    /// Page size for the catalogue scan.
    #[serde(default = "default_page_limit")]
    page_limit: u32,
    /// Safety bound on pages scanned per run.
    #[serde(default = "default_max_pages")]
    max_pages: u32,
    /// Age (seconds) after which a stranded `pending_tier` claim is re-driven (#4).
    #[serde(default = "default_stale_claim_secs")]
    stale_claim_secs: i64,
    /// Queue to enqueue the per-file `tier_move` jobs on.
    #[serde(default)]
    queue: Option<String>,
}

/// Pure due-ness predicate (extracted so it is unit-testable without building a
/// full [`SeaweedFile`]). A file is due for a warm→cold move when it is in the
/// warm tier, its recomputed `cold_move_at` (`reported_at + warm_days`) has
/// elapsed, it is not legally held, and no claim is already in flight.
fn due(
    tier: StorageTier,
    legal_hold: bool,
    pending_tier: Option<StorageTier>,
    reported_at: Option<DateTime<Utc>>,
    warm_days: i64,
    now: DateTime<Utc>,
) -> bool {
    tier == StorageTier::Warm
        && !legal_hold
        && pending_tier.is_none()
        && reported_at
            .map(|r| r + Duration::days(warm_days.max(0)) <= now)
            .unwrap_or(false)
}

/// Pure predicate (#4): whether an in-flight tier-move claim is stale and should
/// be re-driven. A claim is stale when a `pending_tier` is set, the file is not
/// legally held, and the claim is older than `stale_secs` — or has no
/// `pending_since` at all (a pre-`pending_since` claim that can't be aged, treated
/// as stale so it gets cleared rather than stranded forever). Re-driving is safe:
/// `tier_move` is idempotent and CAS-guarded, so it either completes or rolls back.
fn stale_claim(
    pending_tier: Option<StorageTier>,
    pending_since: Option<DateTime<Utc>>,
    legal_hold: bool,
    stale_secs: i64,
    now: DateTime<Utc>,
) -> bool {
    if legal_hold || pending_tier.is_none() {
        return false;
    }
    match pending_since {
        None => true,
        Some(ts) => (now - ts).num_seconds() >= stale_secs,
    }
}

/// Periodic producer: find files due for a warm→cold move and fan out one
/// `tier_move` job per file.
pub struct TierMoveScan {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
}

impl TierMoveScan {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics) -> Self {
        Self { ctx, metrics }
    }

    /// Due = a warm-tier file whose recomputed `cold_move_at` has elapsed, not
    /// legally held, with no claim already in flight.
    fn is_due(file: &SeaweedFile, warm_days: i64, now: DateTime<Utc>) -> bool {
        due(
            file.tier,
            file.legal_hold,
            file.pending_tier,
            file.reported_at,
            warm_days,
            now,
        )
    }

    async fn run_inner(&self, job: Job) -> Result<JobOutcome, WorkerError> {
        let args: TierMoveScanArgs = parse_args(&job)?;
        let api = self.ctx.api()?.clone();
        let queue = args
            .queue
            .clone()
            .unwrap_or_else(|| "lifecycle".to_string());
        let now = Utc::now();

        let mut scanned = 0u64;
        let mut enqueued = 0u64;

        for page in 0..args.max_pages {
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
                // Determine the move to enqueue, if any:
                //  * a fresh warm->cold move when the file is due, or
                //  * a re-drive of a stranded claim to its in-flight target (#4).
                let move_target = if Self::is_due(file, args.warm_days, now) {
                    Some(args.target)
                } else if stale_claim(
                    file.pending_tier,
                    file.pending_since,
                    file.legal_hold,
                    args.stale_claim_secs,
                    now,
                ) {
                    // Re-drive to the claimed target; the idempotent, CAS-guarded
                    // mover either completes the move or rolls the claim back.
                    if file.pending_tier.is_some() {
                        tracing::info!(file = %file.id, "re-driving stale tier-move claim");
                        self.metrics.record(&TelemetryEvent::with_detail(
                            TelemetryOp::TierMove,
                            TelemetryOutcome::Success,
                            "stale_reclaim",
                        ));
                    }
                    file.pending_tier
                } else {
                    None
                };

                let Some(target) = move_target else { continue };

                let job_args = json!({ "file_id": file.id, "target": target });
                // Generous reserve_for: a cold move may run a while; don't let
                // Faktory re-deliver it mid-flight.
                match self
                    .ctx
                    .enqueue(
                        "tier_move",
                        job_args,
                        &queue,
                        Some(StdDuration::from_secs(3600)),
                    )
                    .await
                {
                    Ok(()) => enqueued += 1,
                    Err(e) => tracing::warn!(file = %file.id, "failed to enqueue tier_move: {e}"),
                }
            }

            if (count as u32) < args.page_limit {
                break; // short page = last page
            }
        }

        tracing::info!(scanned, enqueued, target = ?args.target, "tier_move_scan complete");
        self.metrics.record(&TelemetryEvent::with_detail(
            TelemetryOp::TierMove,
            TelemetryOutcome::Success,
            "scan",
        ));
        Ok(JobOutcome::Succeeded)
    }
}

#[async_trait]
impl JobRunner for TierMoveScan {
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

    fn reported(days_ago: i64) -> Option<DateTime<Utc>> {
        Some(Utc::now() - Duration::days(days_ago))
    }

    #[test]
    fn due_warm_file_past_dwell_is_due() {
        assert!(due(
            StorageTier::Warm,
            false,
            None,
            reported(400),
            365,
            Utc::now()
        ));
    }

    #[test]
    fn warm_file_within_dwell_is_not_due() {
        assert!(!due(
            StorageTier::Warm,
            false,
            None,
            reported(10),
            365,
            Utc::now()
        ));
    }

    #[test]
    fn legally_held_file_is_never_due() {
        assert!(!due(
            StorageTier::Warm,
            true,
            None,
            reported(400),
            365,
            Utc::now()
        ));
    }

    #[test]
    fn already_claimed_file_is_not_re_enqueued() {
        assert!(!due(
            StorageTier::Warm,
            false,
            Some(StorageTier::Cold),
            reported(400),
            365,
            Utc::now()
        ));
    }

    #[test]
    fn non_warm_tiers_are_not_due() {
        assert!(!due(
            StorageTier::Hot,
            false,
            None,
            reported(400),
            365,
            Utc::now()
        ));
        assert!(!due(
            StorageTier::Cold,
            false,
            None,
            reported(400),
            365,
            Utc::now()
        ));
    }

    #[test]
    fn file_without_report_out_is_not_due() {
        assert!(!due(StorageTier::Warm, false, None, None, 365, Utc::now()));
    }

    #[test]
    fn tier_label_is_stable() {
        assert_eq!(tier_label(StorageTier::Cold), "cold");
        assert_eq!(tier_label(StorageTier::Warm), "warm");
        assert_eq!(tier_label(StorageTier::Hot), "hot");
    }

    #[test]
    fn no_claim_is_never_stale() {
        assert!(!stale_claim(None, None, false, 3600, Utc::now()));
    }

    #[test]
    fn fresh_claim_is_not_stale() {
        let now = Utc::now();
        assert!(!stale_claim(
            Some(StorageTier::Cold),
            Some(now),
            false,
            3600,
            now
        ));
    }

    #[test]
    fn old_claim_is_stale() {
        let now = Utc::now();
        let since = now - Duration::hours(2);
        assert!(stale_claim(
            Some(StorageTier::Cold),
            Some(since),
            false,
            3600,
            now
        ));
    }

    #[test]
    fn unaged_claim_without_pending_since_is_stale() {
        // pending_tier set but no pending_since (pre-#4 claim) -> treat as stale.
        assert!(stale_claim(
            Some(StorageTier::Cold),
            None,
            false,
            3600,
            Utc::now()
        ));
    }

    #[test]
    fn legally_held_claim_is_never_stale() {
        let now = Utc::now();
        let since = now - Duration::hours(5);
        assert!(!stale_claim(
            Some(StorageTier::Cold),
            Some(since),
            true,
            3600,
            now
        ));
    }
}
