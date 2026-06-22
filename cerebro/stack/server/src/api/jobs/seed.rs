//! Lifecycle schedule seeding (S3-4).
//!
//! Inserts the default periodic [`ScheduleJob`]s that drive the Stage 3 producers —
//! the tier-move scan, retention sweep, purge/reclaim and integrity verify — so a
//! deployment runs the lifecycle on a cadence out of the box, rather than waiting
//! for an operator to hand-create schedules.
//!
//! **Idempotent and operator-respecting.** Each schedule has a *stable* id and is
//! written with `$setOnInsert` under an upsert: the first boot inserts it, later
//! boots leave any existing schedule untouched. So an operator who disables a
//! schedule, retimes it, or edits its args keeps those changes across restarts — the
//! seed only ever *adds what's missing*.
//!
//! Gated by `CEREBRO_SEED_LIFECYCLE_SCHEDULES`; runs only when the in-process
//! scheduler is also enabled (there is no point seeding producers nothing will
//! enqueue).

use chrono::{Duration, Utc};
use mongodb::{
    bson::{doc, to_document},
    options::UpdateOptions,
    Client as MongoClient,
};
use serde_json::json;

use cerebro_model::api::files::retention::RetentionPolicy;
use cerebro_model::api::jobs::model::ScheduleJob;

/// A default schedule to ensure exists.
struct SeedSpec {
    /// Stable id (idempotency key) — never a random UUID, so re-seeding is a no-op.
    id: &'static str,
    kind: &'static str,
    args: serde_json::Value,
    queue: &'static str,
    /// Initial delay from "now" before the first run (staggers the producers).
    offset_secs: i64,
    interval_secs: i64,
    reserve_for_secs: u64,
    retry: i32,
}

/// Seed the default lifecycle schedules into the scheduler's jobs collection.
/// Returns the number of *newly* inserted schedules (existing ones are preserved).
pub async fn seed_lifecycle_schedules(
    mongo: &MongoClient,
    db_name: &str,
    jobs_collection: &str,
) -> anyhow::Result<u32> {
    let coll = mongo
        .database(db_name)
        .collection::<ScheduleJob>(jobs_collection);
    let now = Utc::now();
    let warm_days = RetentionPolicy::from_env().warm_days;
    const DAY: i64 = 86_400;

    // Staggered by a few minutes so the producers don't all fire at once. All daily;
    // operators can retime/disable any of them afterwards and the changes stick.
    let specs = [
        SeedSpec {
            id: "seed:retention_sweep",
            kind: "retention_sweep",
            args: json!({}),
            queue: "maintenance",
            offset_secs: 60,
            interval_secs: DAY,
            reserve_for_secs: 3600,
            retry: 3,
        },
        SeedSpec {
            id: "seed:tier_move_scan",
            kind: "tier_move_scan",
            args: json!({ "target": "Cold", "warm_days": warm_days }),
            queue: "maintenance",
            offset_secs: 120,
            interval_secs: DAY,
            reserve_for_secs: 600,
            retry: 3,
        },
        SeedSpec {
            id: "seed:verify_scan",
            kind: "verify_scan",
            args: json!({ "budget": 500 }),
            queue: "maintenance",
            offset_secs: 180,
            interval_secs: DAY,
            reserve_for_secs: 600,
            retry: 3,
        },
        SeedSpec {
            id: "seed:purge_reclaim",
            kind: "purge_reclaim",
            args: json!({}),
            queue: "maintenance",
            offset_secs: 240,
            interval_secs: DAY,
            reserve_for_secs: 3600,
            retry: 3,
        },
        // Restore recovery (S3-5 #2): a slow hourly safety net that re-drives any
        // archived file stranded mid-restore. Prompt starts come from the
        // restore-request API path; this only backstops dropped poll chains.
        SeedSpec {
            id: "seed:restore_scan",
            kind: "restore_scan",
            args: json!({}),
            queue: "maintenance",
            offset_secs: 300,
            interval_secs: 3600,
            reserve_for_secs: 600,
            retry: 3,
        },
        // Catalogue backup (S4-2): daily mongodump of the control plane to the
        // backup object store. No-op unless the worker has backup settings
        // (CEREBRO_BACKUP_MONGO_URI + CEREBRO_BACKUP_STORE_PATH), so seeding the
        // schedule is safe even before a backup target is configured. Given a
        // generous reservation — a dump + upload of a large catalogue is slower
        // than the lifecycle scans.
        SeedSpec {
            id: "seed:catalogue_backup",
            kind: "catalogue_backup",
            args: json!({}),
            queue: "maintenance",
            offset_secs: 360,
            interval_secs: DAY,
            reserve_for_secs: 6 * 3600,
            retry: 2,
        },
        // Consistency reconcile (S4-3): a weekly, report-first sweep that probes
        // each catalogue entry's backing object and reports dangling references.
        // It never mutates state (S4-5 repairs); orphan reclaim is a separate,
        // operator-gated `reconcile_reclaim` job and is deliberately NOT seeded.
        // Generous reservation — it HEAD-probes the whole estate.
        SeedSpec {
            id: "seed:reconcile_scan",
            kind: "reconcile_scan",
            args: json!({}),
            queue: "maintenance",
            offset_secs: 420,
            interval_secs: 7 * DAY,
            reserve_for_secs: 6 * 3600,
            retry: 2,
        },
        // Verify-repair (S4-5): a weekly pass that consumes the latest reconcile
        // report, confirms each dangling reference against live state (filtering
        // transients and files archived since the scan), and escalates confirmed
        // data loss. Staggered well after reconcile_scan so it reads a fresh report.
        // Mutates nothing, so it is safe to run unattended.
        SeedSpec {
            id: "seed:verify_repair",
            kind: "verify_repair",
            args: json!({}),
            queue: "maintenance",
            offset_secs: 1800,
            interval_secs: 7 * DAY,
            reserve_for_secs: 3 * 3600,
            retry: 2,
        },
        // Archival local-copy reclaim (S4-6 H3): a weekly pass that deletes the
        // redundant local copy of archived files past the grace window, gated on the
        // cold copy being confirmed present. Active only when a cold store is
        // configured; otherwise a no-op. Bounded per run.
        SeedSpec {
            id: "seed:archive_reclaim",
            kind: "archive_reclaim",
            args: json!({}),
            queue: "maintenance",
            offset_secs: 2400,
            interval_secs: 7 * DAY,
            reserve_for_secs: 3 * 3600,
            retry: 2,
        },
    ];

    let total = specs.len() as u32;
    let mut seeded = 0u32;
    for spec in specs {
        let mut job = ScheduleJob::new(
            spec.kind.to_string(),
            spec.args,
            Some(spec.queue.to_string()),
            now + Duration::seconds(spec.offset_secs),
            Some(spec.interval_secs),
            Some(spec.retry),
            Some(spec.reserve_for_secs),
            Some(true),
        );
        job.id = spec.id.to_string();

        // $setOnInsert with the full document minus `id` (the filter supplies it on
        // insert); on a match nothing is written, preserving operator edits.
        let mut bson = to_document(&job)?;
        bson.remove("id");

        let res = coll
            .update_one(doc! { "id": spec.id }, doc! { "$setOnInsert": bson })
            .with_options(UpdateOptions::builder().upsert(true).build())
            .await?;

        if res.upserted_id.is_some() {
            seeded += 1;
            log::info!(
                "Seeded lifecycle schedule '{}' (kind={}, every {}s on queue '{}')",
                spec.id,
                spec.kind,
                spec.interval_secs,
                spec.queue
            );
        }
    }

    log::info!(
        "Lifecycle schedule seeding complete: {} new, {} already present (preserved)",
        seeded,
        total - seeded
    );
    Ok(seeded)
}
