//! Job runners (Stage 3).
//!
//! Each Faktory job `kind` maps to a [`faktory::JobRunner`]. The full lifecycle
//! taxonomy is now implemented:
//! * [`ping::Ping`] — liveness smoke test (S3-1);
//! * [`tier_move`] — `tier_move` + `tier_move_scan` (S3-2a);
//! * [`retention_sweep`] — `retention_sweep` + `purge_reclaim` (S3-2b);
//! * [`restore_drive`] — archival restore state machine (S3-3b);
//! * [`verify`] — `verify_file` + `verify_scan` integrity checks (S3-3a).
//!
//! Every runner records a `started` then a terminal outcome on the worker metrics.

pub mod ping;
pub mod restore_drive;
pub mod retention_sweep;
pub mod tier_move;
pub mod verify;

pub use ping::Ping;
pub use restore_drive::RestoreDrive;
pub use retention_sweep::{PurgeReclaim, RetentionSweep};
pub use tier_move::{TierMove, TierMoveScan};
pub use verify::{VerifyFile, VerifyScan};

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
