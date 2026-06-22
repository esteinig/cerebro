//! Structured telemetry taxonomy for file-lifecycle operations (S2-14).
//!
//! This is **pure data**: it carries no metrics backend. The server maps a
//! [`TelemetryEvent`] onto Prometheus counters (`cerebro-server`'s `telemetry`
//! module), and Stage 3 lifecycle workers can emit the same taxonomy so that
//! operational signals are consistent across the API and the background movers.
//!
//! Label dimensions are deliberately **low-cardinality** ([`TelemetryOp`],
//! [`TelemetryOutcome`], and an optional bounded `detail` such as a tier or
//! retention class) so they are safe to use as Prometheus labels.

use serde::{Deserialize, Serialize};

/// The lifecycle operation a telemetry event describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TelemetryOp {
    /// A file was registered/captured into the platform.
    Upload,
    /// A two-phase storage-tier move (claim → commit).
    TierMove,
    /// A report-out that anchored retention and scheduled placement.
    ReportOut,
    /// A non-destructive expiry (quarantine).
    Expire,
    /// A gated hard purge of quarantined files.
    Purge,
    /// A restore state-machine transition.
    Restore,
    /// An integrity verification pass.
    Verify,
    /// A file deletion.
    Delete,
    /// An audit chain append.
    Audit,
}

impl TelemetryOp {
    /// Stable lower_snake_case label value.
    pub fn as_str(&self) -> &'static str {
        match self {
            TelemetryOp::Upload => "upload",
            TelemetryOp::TierMove => "tier_move",
            TelemetryOp::ReportOut => "report_out",
            TelemetryOp::Expire => "expire",
            TelemetryOp::Purge => "purge",
            TelemetryOp::Restore => "restore",
            TelemetryOp::Verify => "verify",
            TelemetryOp::Delete => "delete",
            TelemetryOp::Audit => "audit",
        }
    }
}

/// The outcome of a lifecycle operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TelemetryOutcome {
    /// The operation completed and changed state as intended.
    Success,
    /// The operation failed unexpectedly (an error condition).
    Failure,
    /// The operation was refused by a guard (precondition / compare-and-set /
    /// protection). Expected and benign — distinguished from `Failure` so alerts
    /// can ignore routine rejections.
    Rejected,
}

impl TelemetryOutcome {
    /// Stable lower_snake_case label value.
    pub fn as_str(&self) -> &'static str {
        match self {
            TelemetryOutcome::Success => "success",
            TelemetryOutcome::Failure => "failure",
            TelemetryOutcome::Rejected => "rejected",
        }
    }
}

/// A single structured telemetry event.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TelemetryEvent {
    pub op: TelemetryOp,
    pub outcome: TelemetryOutcome,
    /// Optional bounded secondary dimension (e.g. a tier `"warm"`, a retention
    /// class, or a restore state). Keep low-cardinality — it becomes a metric
    /// label.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

impl TelemetryEvent {
    /// An event with no secondary dimension.
    pub fn new(op: TelemetryOp, outcome: TelemetryOutcome) -> Self {
        Self {
            op,
            outcome,
            detail: None,
        }
    }

    /// An event with a bounded `detail` dimension.
    pub fn with_detail(
        op: TelemetryOp,
        outcome: TelemetryOutcome,
        detail: impl Into<String>,
    ) -> Self {
        Self {
            op,
            outcome,
            detail: Some(detail.into()),
        }
    }

    /// Convenience: a `Success` event.
    pub fn success(op: TelemetryOp) -> Self {
        Self::new(op, TelemetryOutcome::Success)
    }

    /// Convenience: a `Rejected` (guard-refused) event.
    pub fn rejected(op: TelemetryOp) -> Self {
        Self::new(op, TelemetryOutcome::Rejected)
    }

    /// Convenience: a `Failure` event.
    pub fn failure(op: TelemetryOp) -> Self {
        Self::new(op, TelemetryOutcome::Failure)
    }

    pub fn op_str(&self) -> &'static str {
        self.op.as_str()
    }

    pub fn outcome_str(&self) -> &'static str {
        self.outcome.as_str()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn label_values_are_stable_snake_case() {
        assert_eq!(TelemetryOp::TierMove.as_str(), "tier_move");
        assert_eq!(TelemetryOp::ReportOut.as_str(), "report_out");
        assert_eq!(TelemetryOutcome::Rejected.as_str(), "rejected");
    }

    #[test]
    fn constructors_set_fields() {
        let e =
            TelemetryEvent::with_detail(TelemetryOp::TierMove, TelemetryOutcome::Success, "warm");
        assert_eq!(e.op_str(), "tier_move");
        assert_eq!(e.outcome_str(), "success");
        assert_eq!(e.detail.as_deref(), Some("warm"));

        let r = TelemetryEvent::rejected(TelemetryOp::Restore);
        assert_eq!(r.outcome, TelemetryOutcome::Rejected);
        assert!(r.detail.is_none());
    }
}
