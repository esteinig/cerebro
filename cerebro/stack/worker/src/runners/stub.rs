//! `LifecycleStub` — explicit placeholder for the lifecycle kinds (S3-1).
//!
//! Registered under each lifecycle kind so the worker compiles, connects and
//! exposes the full job taxonomy now, while the real behaviour lands in later
//! packages. It deliberately **fails** (rather than silently succeeding) so a
//! forgotten implementation can never masquerade as a completed lifecycle action.
//!
//! No lifecycle schedules are seeded until the corresponding runner is real
//! (S3-2 / S3-3), so these stubs do not produce retry storms in S3-1 — they only
//! fire if a kind is enqueued by hand.

use std::io;
use std::sync::Arc;

use async_trait::async_trait;
use faktory::{Job, JobRunner};

use crate::context::WorkerContext;
use crate::error::WorkerError;
use crate::telemetry::{JobOutcome, Metrics};

pub struct LifecycleStub {
    kind: &'static str,
    /// Held now so the real runners (S3-2/S3-3) can be dropped in with the same
    /// constructor shape.
    #[allow(dead_code)]
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
}

impl LifecycleStub {
    pub fn new(kind: &'static str, ctx: Arc<WorkerContext>, metrics: Metrics) -> Self {
        Self { kind, ctx, metrics }
    }
}

#[async_trait]
impl JobRunner for LifecycleStub {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {
        self.metrics.record_job(&job.kind().to_string(), JobOutcome::Started);
        tracing::warn!(kind = self.kind, jid = %job.id(), "lifecycle job kind not yet implemented (delivered in a later Stage 3 package)");
        self.metrics.record_job(&job.kind().to_string(), JobOutcome::Failed);
        Err(WorkerError::NotImplemented(self.kind).into())
    }
}
