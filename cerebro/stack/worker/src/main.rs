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

mod config;
mod context;
mod error;
mod health;
mod runners;
mod telemetry;

use faktory::Worker;
use tracing::Level;

use crate::config::WorkerConfig;
use crate::context::WorkerContext;
use crate::runners::{
    Ping, PurgeReclaim, RestoreDrive, RetentionSweep, TierMove, TierMoveScan, VerifyFile, VerifyScan,
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

    // Register runners: a liveness ping plus the full lifecycle taxonomy
    // (movement, retention, restore, integrity).
    let mut builder = Worker::builder();
    builder.register("ping", Ping::new(metrics.clone()));
    builder.register("tier_move", TierMove::new(ctx.clone(), metrics.clone()));
    builder.register("tier_move_scan", TierMoveScan::new(ctx.clone(), metrics.clone()));
    builder.register("retention_sweep", RetentionSweep::new(ctx.clone(), metrics.clone()));
    builder.register("purge_reclaim", PurgeReclaim::new(ctx.clone(), metrics.clone()));
    builder.register("restore_drive", RestoreDrive::new(ctx.clone(), metrics.clone()));
    builder.register("verify_file", VerifyFile::new(ctx.clone(), metrics.clone()));
    builder.register("verify_scan", VerifyScan::new(ctx.clone(), metrics.clone()));

    let mut worker = builder.connect().await.expect("connect to faktory");

    let queues: Vec<&str> = config.queues.iter().map(String::as_str).collect();
    tracing::info!(
        kinds = ?[
            "ping", "tier_move", "tier_move_scan", "retention_sweep",
            "purge_reclaim", "restore_drive", "verify_file", "verify_scan"
        ],
        "registered runners; consuming queues"
    );

    if let Err(e) = worker.run(&queues).await {
        tracing::error!(error = %e, "worker run loop exited with error");
        return Err(error::io_err(e));
    }
    Ok(())
}
