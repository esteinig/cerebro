//! cerebro-worker — Faktory consumer for the Cerebro lifecycle (Stage 3).
//!
//! A standalone, independently-scaled process that pulls jobs from Faktory and
//! executes them. It is the **consumer** half of the producer/consumer split: the
//! API server and its `Scheduler` *enqueue* jobs (on demand or periodically); this
//! process *runs* them. Workers are never spawned inside API request handlers.
//!
//! S3-1 establishes the foundation — shared context (API + FS clients), worker
//! telemetry + health endpoints, and the full job taxonomy — with a real `ping`
//! liveness job and explicit stubs for the lifecycle kinds. Movement workers land
//! in S3-2, integrity/restore workers in S3-3.

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
use crate::runners::{LifecycleStub, Ping, PurgeReclaim, RetentionSweep, TierMove, TierMoveScan};
use crate::telemetry::Metrics;

/// Lifecycle job kinds still registered as stubs (real runners land in S3-3).
/// S3-2a (tier_move/tier_move_scan) and S3-2b (retention_sweep/purge_reclaim) are real.
const LIFECYCLE_KINDS: &[&str] = &[
    "verify_scan",      // S3-3a — scheduled integrity verification + repair
    "restore_drive",    // S3-3b — drive the archival restore state machine
];

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

    // Register runners: a real liveness ping, the S3-2a tier mover + scan, the
    // S3-2b retention sweep + purge/reclaim, and the remaining lifecycle taxonomy
    // as stubs.
    let mut builder = Worker::builder();
    builder.register("ping", Ping::new(metrics.clone()));
    builder.register("tier_move", TierMove::new(ctx.clone(), metrics.clone()));
    builder.register("tier_move_scan", TierMoveScan::new(ctx.clone(), metrics.clone()));
    builder.register("retention_sweep", RetentionSweep::new(ctx.clone(), metrics.clone()));
    builder.register("purge_reclaim", PurgeReclaim::new(ctx.clone(), metrics.clone()));
    for kind in LIFECYCLE_KINDS {
        builder.register(*kind, LifecycleStub::new(*kind, ctx.clone(), metrics.clone()));
    }

    let mut worker = builder.connect().await.expect("connect to faktory");

    let queues: Vec<&str> = config.queues.iter().map(String::as_str).collect();
    tracing::info!(
        real = ?["ping", "tier_move", "tier_move_scan", "retention_sweep", "purge_reclaim"],
        stubs = ?LIFECYCLE_KINDS,
        "registered runners; consuming queues"
    );

    if let Err(e) = worker.run(&queues).await {
        tracing::error!(error = %e, "worker run loop exited with error");
        return Err(error::io_err(e));
    }
    Ok(())
}
