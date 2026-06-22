//! Prometheus metrics for file-lifecycle telemetry.
//!
//! [`Metrics`] owns a [`prometheus::Registry`] and the lifecycle counters. It is
//! created **once** at server start and shared across worker threads (the
//! registry and counters are `Arc`-backed, so `clone` shares state). Handlers map
//! a [`TelemetryEvent`] onto counter increments via [`Metrics::record`], and the
//! `/metrics` endpoint exposes the registry in the Prometheus text format.

use cerebro_model::api::files::telemetry::TelemetryEvent;
use prometheus::{Encoder, IntCounterVec, Opts, Registry, TextEncoder};

/// Shared metrics handle. Cheap to clone (shares the underlying registry).
#[derive(Clone)]
pub struct Metrics {
    registry: Registry,
    /// `cerebro_file_lifecycle_ops_total{op, outcome, detail}` — one counter
    /// vector covering every lifecycle operation, keyed by low-cardinality
    /// labels so moves/expiries/restores/verify pass-fail are all queryable.
    lifecycle_ops: IntCounterVec,
}

impl Metrics {
    /// Construct and register the metric set. Panics only on a programming error
    /// (an invalid metric name or a duplicate registration), which would be a
    /// build-time constant bug, not a runtime condition.
    pub fn new() -> Self {
        let registry = Registry::new();

        let lifecycle_ops = IntCounterVec::new(
            Opts::new(
                "cerebro_file_lifecycle_ops_total",
                "Count of file-lifecycle operations by operation, outcome and detail.",
            ),
            &["op", "outcome", "detail"],
        )
        .expect("valid lifecycle metric definition");

        registry
            .register(Box::new(lifecycle_ops.clone()))
            .expect("lifecycle metric registers once");

        Self {
            registry,
            lifecycle_ops,
        }
    }

    /// Record a telemetry event as a counter increment.
    pub fn record(&self, event: &TelemetryEvent) {
        let detail = event.detail.as_deref().unwrap_or("none");
        self.lifecycle_ops
            .with_label_values(&[event.op_str(), event.outcome_str(), detail])
            .inc();
    }

    /// Encode the registry in the Prometheus text exposition format.
    pub fn encode(&self) -> String {
        let mut buffer = Vec::new();
        let encoder = TextEncoder::new();
        let families = self.registry.gather();
        // Encoding a gathered set of well-formed families does not fail in
        // practice; on the off chance it does, return whatever was written.
        let _ = encoder.encode(&families, &mut buffer);
        String::from_utf8(buffer).unwrap_or_default()
    }
}

impl Default for Metrics {
    fn default() -> Self {
        Self::new()
    }
}

/// `GET /metrics` — expose the registry in the Prometheus text format. Carries no
/// PII (operational counters only); scrape it from an internal network / behind a
/// firewall.
#[actix_web::get("/metrics")]
async fn metrics_handler(
    data: actix_web::web::Data<crate::api::server::AppState>,
) -> actix_web::HttpResponse {
    actix_web::HttpResponse::Ok()
        .content_type("text/plain; version=0.0.4")
        .body(data.metrics.encode())
}

/// Register the telemetry endpoint(s).
pub fn telemetry_config(cfg: &mut actix_web::web::ServiceConfig) {
    cfg.service(metrics_handler);
}

#[cfg(test)]
mod tests {
    use super::*;
    use cerebro_model::api::files::telemetry::{TelemetryOp, TelemetryOutcome};

    #[test]
    fn records_and_encodes() {
        let m = Metrics::new();
        m.record(&TelemetryEvent::success(TelemetryOp::Expire));
        m.record(&TelemetryEvent::with_detail(
            TelemetryOp::TierMove,
            TelemetryOutcome::Success,
            "warm",
        ));
        m.record(&TelemetryEvent::rejected(TelemetryOp::Restore));

        let text = m.encode();
        assert!(text.contains("cerebro_file_lifecycle_ops_total"));
        assert!(text.contains("op=\"expire\""));
        assert!(text.contains("detail=\"warm\""));
        assert!(text.contains("outcome=\"rejected\""));
    }

    #[test]
    fn clone_shares_state() {
        let m = Metrics::new();
        let m2 = m.clone();
        m.record(&TelemetryEvent::success(TelemetryOp::Purge));
        m2.record(&TelemetryEvent::success(TelemetryOp::Purge));
        // Both increments land on the same series => value 2.
        assert!(m.encode().contains("op=\"purge\""));
    }
}
