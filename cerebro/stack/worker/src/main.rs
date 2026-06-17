//! cerebro-worker — Faktory consumer for the Cerebro lifecycle (Stage 3).
//!
//! A standalone, independently-scaled process that pulls jobs from Faktory and
//! executes them. It is the **consumer** half of the producer/consumer split: the
//! API server and its `Scheduler` *enqueue* jobs (on demand or periodically); this
//! process *runs* them. Workers are never spawned inside API request handlers.
//!
//! S3-1 established the foundation — shared context (API + FS clients), worker
//! telemetry + health endpoints, and the job taxonomy. The lifecycle runners are
//! now all real: tier mover (S3-2a), retention sweep + purge/reclaim (S3-2b),
//! restore executor (S3-3b), and integrity verify (S3-3a).

mod archive;
mod backup;
mod config;
mod context;
mod error;
mod health;
mod reconcile;
mod runners;
mod telemetry;

use faktory::Worker;
use tracing::Level;

use crate::config::WorkerConfig;
use crate::context::WorkerContext;
use crate::runners::{
    CatalogueBackup, Ping, PurgeReclaim, ReconcileReclaim, ReconcileScan, RestoreDrive, RestoreScan, RetentionSweep, TierMove, TierMoveScan, VerifyFile, VerifyRepair, VerifyScan, ArchiveReclaim,
};
use crate::telemetry::Metrics;

fn init_tracing() {
    tracing_subscriber::fmt()
        .with_env_filter(std::env::var("RUST_LOG").unwrap_or_else(|_| "info".into()))
        .with_level(true)
        .with_target(false)
        .with_max_level(Level::INFO)
        .init();
}

#[tokio::main]
async fn main() -> std::io::Result<()> {
    init_tracing();

    let config = WorkerConfig::from_env();
    tracing::info!(queues = ?config.queues, metrics = %config.metrics_addr, "cerebro-worker starting");
    match std::env::var("FAKTORY_URL") {
        Ok(url) => tracing::info!(faktory = %url, "FAKTORY_URL set"),
        Err(_) => tracing::warn!("FAKTORY_URL not set; using faktory client defaults"),
    }

    // Worker telemetry + health/metrics endpoint.
    let metrics = Metrics::new();
    health::spawn(metrics.clone(), config.metrics_addr.clone())?;

    // Build the shared context (blocking reqwest clients) off the async executor.
    let ctx = {
        let cfg = config.clone();
        tokio::task::spawn_blocking(move || WorkerContext::build(cfg))
            .await
            .map_err(error::io_err)?
    };

    // Authenticate as the service Bot (S3-5 #5). If bot credentials are configured
    // this logs in and rebuilds the clients with a long-lived, role-scoped token; a
    // failure is logged but non-fatal (lifecycle runners will error until auth
    // succeeds — e.g. on a later periodic re-login or restart).
    if let Err(e) = ctx.login().await {
        tracing::error!(error = %e, "initial bot login failed; lifecycle runners will error until re-login succeeds");
    }

    // Periodic re-login keeps the token fresh well within its TTL and rebuilds the
    // clients in place (decision #3). No-op when bot credentials aren't set.
    if config.bot_email.is_some() {
        let relogin_ctx = ctx.clone();
        let interval = config.relogin_interval_secs.max(60);
        tokio::spawn(async move {
            let mut tick = tokio::time::interval(std::time::Duration::from_secs(interval));
            tick.tick().await; // consume the immediate first tick (already logged in)
            loop {
                tick.tick().await;
                if let Err(e) = relogin_ctx.relogin().await {
                    tracing::warn!(error = %e, "periodic bot re-login failed; will retry next interval");
                }
            }
        });
    }

    // Register runners: a liveness ping plus the full lifecycle taxonomy
    // (movement, retention, restore, integrity).
        // Register runners: a liveness ping plus the full lifecycle taxonomy
    // (movement, retention, restore, integrity).
    let mut worker = Worker::builder()
        .register("ping", Ping::new(metrics.clone()))
        .register("tier_move", TierMove::new(ctx.clone(), metrics.clone(), config.archive.clone()))
        .register("tier_move_scan", TierMoveScan::new(ctx.clone(), metrics.clone()))
        .register("retention_sweep", RetentionSweep::new(ctx.clone(), metrics.clone()))
        .register("purge_reclaim", PurgeReclaim::new(ctx.clone(), metrics.clone()))
        .register("restore_drive", RestoreDrive::new(ctx.clone(), metrics.clone(), config.archive.clone()))
        .register("restore_scan", RestoreScan::new(ctx.clone(), metrics.clone()))
        .register("verify_file", VerifyFile::new(ctx.clone(), metrics.clone(), config.archive.clone()))
        .register("verify_scan", VerifyScan::new(ctx.clone(), metrics.clone()))
        .register("catalogue_backup", CatalogueBackup::new(ctx.clone(), metrics.clone(), config.backup.clone()))
        .register("reconcile_scan", ReconcileScan::new(ctx.clone(), metrics.clone(), config.backup.as_ref().map(|b| b.store_root.clone())))
        .register("reconcile_reclaim", ReconcileReclaim::new(ctx.clone(), metrics.clone()))
        .register("verify_repair", VerifyRepair::new(ctx.clone(), metrics.clone(), config.backup.as_ref().map(|b| b.store_root.clone())))
        .register("archive_reclaim", ArchiveReclaim::new(ctx.clone(), metrics.clone(), config.archive.clone()))
        .connect()
        .await
        .expect("connect to faktory");

    let queues: Vec<&str> = config.queues.iter().map(String::as_str).collect();
    tracing::info!(
        kinds = ?[
            "ping", "tier_move", "tier_move_scan", "retention_sweep",
            "purge_reclaim", "restore_drive", "restore_scan", "verify_file", "verify_scan"
        ],
        "registered runners; consuming queues"
    );

    if let Err(e) = worker.run(&queues).await {
        tracing::error!(error = %e, "worker run loop exited with error");
        return Err(error::io_err(e));
    }
    Ok(())
}
