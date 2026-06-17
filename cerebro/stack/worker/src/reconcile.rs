//! Consistency reconcile engine (S4-3).
//!
//! Detects divergence between the catalogue (MongoDB, the authoritative record of
//! what *should* exist) and the object store (SeaweedFS, what *does* exist):
//!
//! * **Dangling references** — a catalogue entry whose backing object is missing.
//!   This is latent data loss; the scan reports it (and S4-5 repairs from a
//!   replica/backup). Safe to detect automatically.
//! * **Orphan objects** — a stored object with no catalogue entry. Wasted space,
//!   and deleting one is destructive, so reclaim is operator-gated (D6).
//!
//! Both directions are pure set differences with a **grace window** so the scan
//! never acts on a just-written object whose counterpart is moments behind. The
//! functions here are deterministic and unit-tested; the runner supplies the two
//! sides (catalogue via the API, object presence via `FileSystemClient`).

use std::collections::HashSet;

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

/// A catalogue entry projected to what reconcile needs.
#[derive(Debug, Clone)]
pub struct CatalogueRef {
    /// Catalogue document id.
    pub file_id: String,
    /// Backing object id (SeaweedFS fid).
    pub fid: String,
    /// Registration time, when parseable — used for grace filtering.
    pub created_at: Option<DateTime<Utc>>,
    /// Whether the object has been archived to remote cold (S3 Glacier). Archived
    /// objects are legitimately absent from the local store, so they are never
    /// dangling — retrieval goes through the restore path, not reconcile.
    pub archived: bool,
    /// Current storage tier label (for the report).
    pub tier: String,
}

/// An object found in the store.
#[derive(Debug, Clone)]
pub struct StoreObjectRef {
    /// Store key (fid or filer path) — compared against catalogue fids.
    pub key: String,
    /// Last-modified time, when known — used for grace filtering.
    pub modified_at: Option<DateTime<Utc>>,
}

/// A catalogue entry whose backing object is missing (latent data loss).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct DanglingRef {
    pub file_id: String,
    pub fid: String,
    pub tier: String,
}

/// A stored object with no catalogue entry (reclaim candidate).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct OrphanRef {
    pub key: String,
}

/// The result of a reconcile pass — the report-first output (D6).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReconcileReport {
    pub generated_at: DateTime<Utc>,
    /// Catalogue entries considered this pass.
    pub catalogue_count: usize,
    /// Entries whose object presence was actually probed (archived/errored
    /// entries are excluded).
    pub probed_count: usize,
    /// Whether the store side was independently enumerated. Orphan detection is
    /// only meaningful when this is `true`; until store enumeration lands it is
    /// `false` and `orphans` is empty (the operator-gated reclaim still accepts an
    /// explicit key list).
    pub store_enumerated: bool,
    pub store_count: Option<usize>,
    /// Probes that errored (treated as unknown, never as missing).
    pub probe_errors: usize,
    pub dangling: Vec<DanglingRef>,
    pub orphans: Vec<OrphanRef>,
}

/// True when `ts` is known and within `grace` of `now` — i.e. too recent to act
/// on. An unknown timestamp is *not* grace-protected (returns `false`): the
/// catalogue records objects only after a successful upload, so a missing object
/// for an entry of unknown age is still a real finding.
fn within_grace(ts: Option<DateTime<Utc>>, grace: chrono::Duration, now: DateTime<Utc>) -> bool {
    matches!(ts, Some(t) if now - t < grace)
}

/// Catalogue entries whose backing object is absent from `present_fids`.
///
/// Skips archived entries (their object lives in remote cold, not the local store)
/// and entries newer than `grace` (an in-flight registration race). Whatever
/// remains and is not in `present_fids` is dangling.
pub fn find_dangling(
    catalogue: &[CatalogueRef],
    present_fids: &HashSet<String>,
    grace: chrono::Duration,
    now: DateTime<Utc>,
) -> Vec<DanglingRef> {
    catalogue
        .iter()
        .filter(|c| !c.archived)
        .filter(|c| !within_grace(c.created_at, grace, now))
        .filter(|c| !present_fids.contains(&c.fid))
        .map(|c| DanglingRef { file_id: c.file_id.clone(), fid: c.fid.clone(), tier: c.tier.clone() })
        .collect()
}

/// Stored objects with no catalogue entry.
///
/// Skips objects newer than `grace` (a just-written object whose catalogue entry
/// is moments behind). Whatever remains and is not a known catalogue fid is an
/// orphan reclaim candidate (acted on only via the operator-gated reclaim path).
pub fn find_orphans(
    store: &[StoreObjectRef],
    catalogue_fids: &HashSet<String>,
    grace: chrono::Duration,
    now: DateTime<Utc>,
) -> Vec<OrphanRef> {
    store
        .iter()
        .filter(|o| !within_grace(o.modified_at, grace, now))
        .filter(|o| !catalogue_fids.contains(&o.key))
        .map(|o| OrphanRef { key: o.key.clone() })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::Duration;

    fn cat(file_id: &str, fid: &str, age_days: Option<i64>, archived: bool, now: DateTime<Utc>) -> CatalogueRef {
        CatalogueRef {
            file_id: file_id.into(),
            fid: fid.into(),
            created_at: age_days.map(|d| now - Duration::days(d)),
            archived,
            tier: "Hot".into(),
        }
    }

    #[test]
    fn dangling_flags_missing_old_local_objects_only() {
        let now = Utc::now();
        let catalogue = vec![
            cat("f1", "3,aaa", Some(30), false, now), // present -> not dangling
            cat("f2", "3,bbb", Some(30), false, now), // absent + old + local -> DANGLING
            cat("f3", "3,ccc", Some(30), true, now),  // absent but archived -> skip
            cat("f4", "3,ddd", Some(0), false, now),  // absent but too new (grace) -> skip
            cat("f5", "3,eee", None, false, now),     // absent, unknown age -> DANGLING
        ];
        let present: HashSet<String> = ["3,aaa".to_string()].into_iter().collect();
        let dangling = find_dangling(&catalogue, &present, Duration::days(1), now);
        let ids: Vec<&str> = dangling.iter().map(|d| d.file_id.as_str()).collect();
        assert_eq!(ids, vec!["f2", "f5"]);
    }

    #[test]
    fn orphans_are_store_minus_catalogue_past_grace() {
        let now = Utc::now();
        let store = vec![
            StoreObjectRef { key: "3,aaa".into(), modified_at: Some(now - Duration::days(10)) }, // catalogued
            StoreObjectRef { key: "3,zzz".into(), modified_at: Some(now - Duration::days(10)) }, // ORPHAN
            StoreObjectRef { key: "3,fresh".into(), modified_at: Some(now) },                    // too new -> skip
            StoreObjectRef { key: "3,unknown".into(), modified_at: None },                       // ORPHAN (no grace)
        ];
        let catalogue_fids: HashSet<String> = ["3,aaa".to_string()].into_iter().collect();
        let orphans = find_orphans(&store, &catalogue_fids, Duration::days(1), now);
        let keys: Vec<&str> = orphans.iter().map(|o| o.key.as_str()).collect();
        assert_eq!(keys, vec!["3,zzz", "3,unknown"]);
    }

    #[test]
    fn within_grace_handles_known_and_unknown() {
        let now = Utc::now();
        assert!(within_grace(Some(now), Duration::days(1), now));            // just now -> protected
        assert!(!within_grace(Some(now - Duration::days(2)), Duration::days(1), now)); // old -> not
        assert!(!within_grace(None, Duration::days(1), now));               // unknown -> not protected
    }
}
