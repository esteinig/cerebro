//! Worker-side telemetry.
//!
//! Two metric families, exposed at the worker's own `/metrics`:
//!
//! * `cerebro_worker_jobs_total{kind, outcome}` — every job the worker runs, by
//!   kind and terminal outcome. This is the worker's own health signal.
//! * `cerebro_file_lifecycle_ops_total{op, outcome, detail}` — the **same** family
//!   the server exposes, for lifecycle effects the worker performs that do
//!   **not** flow through an API endpoint (e.g. a `verify` pass/fail, bytes
//!   reclaimed). API-mediated effects are already counted server-side, so they are
//!   not double-counted here.

use cerebro_model::api::files::telemetry::TelemetryEvent;
use prometheus::{Encoder, IntCounterVec, Opts, Registry, TextEncoder};

/// Terminal outcome of a job run.
#[derive(Debug, Clone, Copy)]
pub enum JobOutcome {
    Started,
    Succeeded,
    Failed,
    /// Ran but intentionally did nothing (e.g. a not-yet-implemented stub).
    Skipped,
}

impl JobOutcome {
    fn as_str(&self) -> &'static str {
        match self {
            JobOutcome::Started => "started",
            JobOutcome::Succeeded => "succeeded",
            JobOutcome::Failed => "failed",
            JobOutcome::Skipped => "skipped",
        }
    }
}

/// Shared metrics handle. Cheap to clone (shares the underlying registry).
#[derive(Clone)]
pub struct Metrics {
    registry: Registry,
    jobs: IntCounterVec,
    lifecycle: IntCounterVec,
}

impl Metrics {
    pub fn new() -> Self {
        let registry = Registry::new();

        let jobs = IntCounterVec::new(
            Opts::new(
                "cerebro_worker_jobs_total",
                "Worker jobs run, by kind and outcome.",
            ),
            &["kind", "outcome"],
        )
        .expect("valid worker jobs metric");

        let lifecycle = IntCounterVec::new(
            Opts::new(
                "cerebro_file_lifecycle_ops_total",
                "Worker-performed file-lifecycle operations by op, outcome and detail.",
            ),
            &["op", "outcome", "detail"],
        )
        .expect("valid lifecycle metric");

        registry
            .register(Box::new(jobs.clone()))
            .expect("register jobs");
        registry
            .register(Box::new(lifecycle.clone()))
            .expect("register lifecycle");

        Self {
            registry,
            jobs,
            lifecycle,
        }
    }

    /// Record a job run outcome.
    pub fn record_job(&self, kind: &str, outcome: JobOutcome) {
        self.jobs.with_label_values(&[kind, outcome.as_str()]).inc();
    }

    /// Record a worker-performed lifecycle event (not API-mediated).
    pub fn record(&self, event: &TelemetryEvent) {
        let detail = event.detail.as_deref().unwrap_or("none");
        self.lifecycle
            .with_label_values(&[event.op_str(), event.outcome_str(), detail])
            .inc();
    }

    /// Encode the registry in the Prometheus text exposition format.
    pub fn encode(&self) -> String {
        let mut buffer = Vec::new();
        let encoder = TextEncoder::new();
        let families = self.registry.gather();
        let _ = encoder.encode(&families, &mut buffer);
        String::from_utf8(buffer).unwrap_or_default()
    }
}

impl Default for Metrics {
    fn default() -> Self {
        Self::new()
    }
}
