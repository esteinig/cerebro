use crate::api::{files::model::SeaweedFileId, watchers::model::ProductionWatcher};
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

use super::model::{FileTag, FileType};
use super::retention::{RestoreState, RetentionClass, StorageTier};

/// Dedicated archive/relocate request.
///
/// Repoints a file between local and remote (archival) storage in one
/// compare-and-set operation, kept separate from [`UpdateFileLifecycleSchema`] so
/// the fid repoint — which atomically changes `tier`, `archived`, `fid`, and
/// `archive_key` — is never tangled with routine lifecycle edits. Applied only if
/// the file's current tier equals `expected_tier`, so concurrent movers cannot
/// double-apply.
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct FileRelocateSchema {
    /// CAS guard: apply only if the file's current `tier` matches.
    pub expected_tier: StorageTier,
    /// Tier after relocation (`Cold` when archiving; a local tier when restoring).
    pub target_tier: StorageTier,
    /// Archived flag after relocation: `true` once the bytes live in the remote
    /// store, `false` once restored to local.
    pub archived: bool,
    /// New backing fid. Set on **restore** to the freshly re-uploaded local
    /// object; `None` leaves the existing fid (e.g. archiving in place before the
    /// local copy is reclaimed).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fid: Option<String>,
    /// Remote object key for the archived bytes. Set on **archive**, cleared on
    /// **restore** (send `clear_archive_key`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub archive_key: Option<String>,
    /// Clear `archive_key` (used on restore). Takes precedence over `archive_key`.
    #[serde(default)]
    pub clear_archive_key: bool,
}

#[derive(Deserialize, Serialize, Debug)]
pub struct RegisterFileSchema {
    pub id: String,
    pub date: String,
    pub name: String,
    pub hash: String,
    pub size: u64,
    pub fid: SeaweedFileId,
    pub run_id: Option<String>,
    pub sample_id: Option<String>,
    pub pipeline_id: Option<String>,
    pub description: Option<String>,
    pub ftype: Option<FileType>,
    pub watcher: Option<ProductionWatcher>,
    /// Filer object path when stored via the path-addressed filer.
    #[serde(default)]
    pub path: Option<String>,
    /// Physical storage tier at registration.
    #[serde(default)]
    pub tier: StorageTier,
    /// Retention category assigned to the file.
    #[serde(default)]
    pub retention: RetentionClass,
    /// Absolute expiry computed by the registering client from its retention
    /// policy. `None` means "retain indefinitely".
    #[serde(default)]
    pub retain_until: Option<DateTime<Utc>>,
    /// When set, the file is exempt from expiry.
    #[serde(default)]
    pub legal_hold: bool,
    /// Requested/observed replica count, when known.
    #[serde(default)]
    pub replicas: Option<u32>,
    /// Whether the object's data resides in remote archival storage (Glacier)
    /// and requires a restore before retrieval.
    #[serde(default)]
    pub archived: bool,
    /// When the result for this file's case was reported out (retention anchor).
    #[serde(default)]
    pub reported_at: Option<DateTime<Utc>>,
}

#[derive(Deserialize, Serialize, Debug)]
pub struct UpdateFileTagsSchema {
    pub ids: Vec<String>,
    pub tags: Vec<FileTag>,
}

/// Partial update of a file's lifecycle fields.
///
/// Every field is optional: one left unset (`None`/`false`) is **unchanged**.
/// This backs the `PATCH /files/{id}/lifecycle` endpoint and lets the report-out
/// (FS-7/FS-9), restore (FS-4), and integrity (FS-6) flows persist their computed
/// state without re-registering the file. Placement (`tier`) and compliance
/// (`retain_until`, `legal_hold`) remain independent axes.
#[derive(Deserialize, Serialize, Debug, Default, Clone)]
pub struct UpdateFileLifecycleSchema {
    /// Move the file to this storage tier (placement).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tier: Option<StorageTier>,
    /// Set the retention expiry to this timestamp (compliance).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub retain_until: Option<DateTime<Utc>>,
    /// Clear the retention expiry (retain indefinitely). Takes precedence over
    /// `retain_until` when both are supplied.
    #[serde(default)]
    pub clear_retain_until: bool,
    /// Record the report-out timestamp (the retention anchor).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reported_at: Option<DateTime<Utc>>,
    /// Place (`true`) or release (`false`) a legal hold.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub legal_hold: Option<bool>,
    /// Mark the object archived on Glacier (`true`) or not (`false`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub archived: Option<bool>,
    /// Update the archival restore state.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub restore_state: Option<RestoreState>,
    /// Update the observed replica count.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub replicas: Option<u32>,
    /// Claim an in-flight tier move to this target. Set during the move;
    /// cleared automatically when `tier` is committed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pending_tier: Option<StorageTier>,
    /// Clear the pending tier claim (e.g. aborting a move). Takes precedence over
    /// `pending_tier`.
    #[serde(default)]
    pub clear_pending_tier: bool,
    /// Compare-and-set precondition: the update applies only if the file's
    /// current `tier` equals this. Makes mover updates idempotent under
    /// at-least-once delivery. Not itself a mutation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_tier: Option<StorageTier>,
    /// Stamp the claim time for an in-flight tier move. The mover sets
    /// this alongside `pending_tier` so a crashed move can be aged out by the scan.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pending_since: Option<DateTime<Utc>>,
    /// Clear `pending_since` (paired with `clear_pending_tier` on commit/rollback).
    #[serde(default)]
    pub clear_pending_since: bool,
    /// Stamp a successful integrity verification. The verify runner sets
    /// this so the scan can order by least-recently-verified for full coverage.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub verified_at: Option<DateTime<Utc>>,
}

impl UpdateFileLifecycleSchema {
    /// True when the update would change nothing — used by the endpoint to reject
    /// no-op requests. `expected_tier` is a precondition, not a mutation, so it
    /// does not count.
    pub fn is_empty(&self) -> bool {
        self.tier.is_none()
            && self.retain_until.is_none()
            && !self.clear_retain_until
            && self.reported_at.is_none()
            && self.legal_hold.is_none()
            && self.archived.is_none()
            && self.restore_state.is_none()
            && self.replicas.is_none()
            && self.pending_tier.is_none()
            && !self.clear_pending_tier
            && self.pending_since.is_none()
            && !self.clear_pending_since
            && self.verified_at.is_none()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lifecycle_update_partial_deserialize() {
        // Only two fields supplied; the rest must remain unset.
        let json = r#"{"tier":"Cold","legal_hold":true}"#;
        let schema: UpdateFileLifecycleSchema = serde_json::from_str(json).unwrap();
        assert_eq!(schema.tier, Some(StorageTier::Cold));
        assert_eq!(schema.legal_hold, Some(true));
        assert!(schema.retain_until.is_none());
        assert!(schema.reported_at.is_none());
        assert!(!schema.clear_retain_until);
        assert!(!schema.is_empty());
    }

    #[test]
    fn lifecycle_update_empty_is_noop() {
        let schema: UpdateFileLifecycleSchema = serde_json::from_str("{}").unwrap();
        assert!(schema.is_empty());
    }

    #[test]
    fn clear_retain_until_is_not_empty() {
        let schema = UpdateFileLifecycleSchema {
            clear_retain_until: true,
            ..Default::default()
        };
        assert!(!schema.is_empty());
    }
}

/// Trigger the report-out lifecycle for a run (optionally a single sample).
///
/// Backs the server-authoritative `POST /files/report-out`: the server applies
/// its configured retention policy to every matching file, anchoring retention at
/// `reported_at` and moving the data off the hot tier (FS-7/FS-9).
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct ReportOutSchema {
    pub run_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    /// Report-out timestamp; defaults to now server-side when omitted.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reported_at: Option<DateTime<Utc>>,
}

/// Trigger a non-destructive expiry sweep: quarantine files whose
/// retention has lapsed (and which are not under legal hold), optionally scoped
/// to a run/sample. `dry_run` previews the eligible set without changing state.
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct ExpireSchema {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    #[serde(default)]
    pub dry_run: bool,
}

/// Trigger a gated hard purge: permanently remove quarantined files that
/// have cleared the quarantine grace window (and are not under legal hold).
/// `dry_run` previews the eligible set without deleting anything.
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct PurgeSchema {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    #[serde(default)]
    pub dry_run: bool,
}

/// Drive a restore state-machine transition.
///
/// The restore worker advances `target` (e.g. `Requested` → `InProgress` →
/// `Restored`); the server validates the transition and stamps the relevant
/// timestamps.
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct RestoreTransitionSchema {
    pub target: RestoreState,
}
