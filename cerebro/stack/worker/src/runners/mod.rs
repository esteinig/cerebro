//! Job runners (S3-1).
//!
//! Each Faktory job `kind` maps to a [`faktory::JobRunner`]. S3-1 ships:
//! * [`ping::Ping`] — a real liveness job (end-to-end smoke test), and
//! * [`stub::LifecycleStub`] — registered under every lifecycle kind, surfacing an
//!   explicit "not yet implemented" error so the seam is honest. The real bodies
//!   land in S3-2 (movement) and S3-3 (integrity/restore).
//!
//! Every runner records a `started` then a terminal outcome on the worker metrics.

pub mod ping;
pub mod retention_sweep;
pub mod stub;
pub mod tier_move;

pub use ping::Ping;
pub use retention_sweep::{PurgeReclaim, RetentionSweep};
pub use stub::LifecycleStub;
pub use tier_move::{TierMove, TierMoveScan};

use faktory::Job;
use serde::Deserialize;

use crate::error::WorkerError;

/// Pull the first JSON object from `job.args()` and deserialize into `T`.
pub fn parse_args<T: for<'de> Deserialize<'de>>(job: &Job) -> Result<T, WorkerError> {
    let first = job
        .args()
        .first()
        .ok_or_else(|| WorkerError::InvalidArgs("missing job args (expected a single JSON object)".into()))?;
    serde_json::from_value::<T>(first.clone())
        .map_err(|e| WorkerError::InvalidArgs(e.to_string()))
}
