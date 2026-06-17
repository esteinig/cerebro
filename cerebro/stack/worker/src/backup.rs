//! Catalogue backup engine (S4-2).
//!
//! The control plane — MongoDB holding the file catalogue, lifecycle state,
//! provenance, and the tamper-evident audit chain — is irreplaceable. Lose it and
//! every stored object is orphaned and the audit history is gone. This module
//! provides the pieces the `catalogue_backup` runner composes:
//!
//! * [`ObjectStore`] — a swappable backup destination (D1). The default is
//!   [`FilesystemObjectStore`] (a mounted backup volume / NFS); an S3-compatible
//!   backend lands in S4-4 behind this same trait.
//! * [`BackupManifest`] — a verifiable record of each backup (checksum, sizes,
//!   and whether the audit chain verified intact at backup time).
//! * [`select_for_deletion`] — retention selection (keep-last-N with an age
//!   floor), kept pure for unit testing.

use std::path::{Path, PathBuf};

use chrono::{DateTime, NaiveDateTime, Utc};
use serde::{Deserialize, Serialize};

/// Manifest schema version, bumped if the on-disk format changes so a restore can
/// reject a manifest it does not understand.
pub const MANIFEST_FORMAT_VERSION: u32 = 1;

/// Timestamp format used for backup ids and object key prefixes (sortable, no
/// separators that are awkward in object keys).
pub const BACKUP_ID_FORMAT: &str = "%Y%m%dT%H%M%SZ";

/// A pluggable backup destination (D1).
///
/// Keys are forward-slash-joined paths under the store's root/prefix. The backup
/// engine is written entirely against this trait, so the target is swappable
/// (filesystem now; AWS S3 / MinIO / SeaweedFS-S3 in S4-4) without touching the
/// backup or retention logic.
pub trait ObjectStore: Send + Sync {
    /// Store `bytes` at `key`, overwriting any existing object.
    fn put(&self, key: &str, bytes: &[u8]) -> anyhow::Result<()>;
    /// Fetch the object stored at `key`.
    fn get(&self, key: &str) -> anyhow::Result<Vec<u8>>;
    /// List every key that begins with `prefix`.
    fn list(&self, prefix: &str) -> anyhow::Result<Vec<String>>;
    /// Delete the object at `key`. Succeeds (no-op) if the key is already absent.
    fn delete(&self, key: &str) -> anyhow::Result<()>;
}

/// Filesystem-backed object store: objects are files under `root`.
///
/// The default backend, and a legitimate independent target as long as `root` is
/// on a disk physically separate from the SeaweedFS data (or, for host-failure
/// durability, an NFS / remote mount). See `docs/catalogue-backup.md`.
pub struct FilesystemObjectStore {
    root: PathBuf,
}

impl FilesystemObjectStore {
    pub fn new(root: impl Into<PathBuf>) -> Self {
        Self { root: root.into() }
    }

    /// Map an object key onto a path under `root`, dropping empty and `.`/`..`
    /// segments so a key can never escape the store root.
    fn path_for(&self, key: &str) -> PathBuf {
        let mut p = self.root.clone();
        for seg in key.split('/').filter(|s| !s.is_empty() && *s != "." && *s != "..") {
            p.push(seg);
        }
        p
    }
}

impl ObjectStore for FilesystemObjectStore {
    fn put(&self, key: &str, bytes: &[u8]) -> anyhow::Result<()> {
        let path = self.path_for(key);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        std::fs::write(&path, bytes)?;
        Ok(())
    }

    fn get(&self, key: &str) -> anyhow::Result<Vec<u8>> {
        Ok(std::fs::read(self.path_for(key))?)
    }

    fn list(&self, prefix: &str) -> anyhow::Result<Vec<String>> {
        fn walk(dir: &Path, root: &Path, out: &mut Vec<String>) -> anyhow::Result<()> {
            if !dir.exists() {
                return Ok(());
            }
            for entry in std::fs::read_dir(dir)? {
                let path = entry?.path();
                if path.is_dir() {
                    walk(&path, root, out)?;
                } else if let Ok(rel) = path.strip_prefix(root) {
                    out.push(rel.to_string_lossy().replace(std::path::MAIN_SEPARATOR, "/"));
                }
            }
            Ok(())
        }
        let mut keys = Vec::new();
        walk(&self.root, &self.root, &mut keys)?;
        keys.retain(|k| k.starts_with(prefix));
        keys.sort();
        Ok(keys)
    }

    fn delete(&self, key: &str) -> anyhow::Result<()> {
        match std::fs::remove_file(self.path_for(key)) {
            Ok(()) => Ok(()),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(e) => Err(e.into()),
        }
    }
}

/// A verifiable record written alongside every backup archive.
///
/// The restore path reads this first: it checks `format_version`, verifies the
/// archive against `archive_blake3`, and — crucially — re-verifies the audit
/// chain after loading and compares it against `audit_chain_verified`, so a
/// backup cannot silently launder a broken chain.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct BackupManifest {
    pub format_version: u32,
    /// Sortable id and key prefix, e.g. `20260617T031500Z`.
    pub backup_id: String,
    pub created_at: DateTime<Utc>,
    pub mongo_db: String,
    /// Object key of the gzip archive within the store.
    pub archive_key: String,
    pub archive_bytes: u64,
    /// BLAKE3 hex digest of the archive (matches `cerebro-fs` hashing).
    pub archive_blake3: String,
    /// Whether the tamper-evident audit chain verified intact at backup time.
    /// `None` when the check could not run (e.g. the API was unreachable); the
    /// restore surfaces this rather than treating it as verified.
    pub audit_chain_verified: Option<bool>,
    /// Number of audit events observed during the pre-backup chain check.
    pub audit_event_count: usize,
    /// `cerebro` version that produced the backup (best-effort provenance).
    pub cerebro_version: String,
}

/// A backup already present in the store, identified by its id and creation time.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BackupRef {
    pub backup_id: String,
    pub created_at: DateTime<Utc>,
}

/// Recover [`BackupRef`]s from the manifest keys under `prefix`.
///
/// Each backup lives at `<prefix>/<backup_id>/manifest.json`; the id is itself a
/// timestamp ([`BACKUP_ID_FORMAT`]), so the creation time is parsed from the id
/// without fetching each manifest. Keys that don't match the layout are ignored.
pub fn parse_backup_refs(prefix: &str, keys: &[String]) -> Vec<BackupRef> {
    let stem = format!("{}/", prefix.trim_end_matches('/'));
    let mut refs = Vec::new();
    for k in keys {
        let Some(rest) = k.strip_prefix(&stem) else { continue };
        let Some(id) = rest.strip_suffix("/manifest.json") else { continue };
        if let Ok(naive) = NaiveDateTime::parse_from_str(id, BACKUP_ID_FORMAT) {
            refs.push(BackupRef {
                backup_id: id.to_string(),
                created_at: DateTime::<Utc>::from_naive_utc_and_offset(naive, Utc),
            });
        }
    }
    refs
}

/// Select which existing backups to delete under a keep-last-N policy with an age
/// floor (S4-2 / D5).
///
/// The `keep_last` most-recent backups are always retained. Older ones are
/// eligible for deletion, but only once they are at least `min_age` old — so a
/// burst of manual backups on one day cannot immediately evict the daily history
/// that provides point-in-time coverage. Pure and deterministic for testing; the
/// runner deletes whatever this returns.
pub fn select_for_deletion(
    mut existing: Vec<BackupRef>,
    keep_last: usize,
    min_age: chrono::Duration,
    now: DateTime<Utc>,
) -> Vec<BackupRef> {
    existing.sort_by(|a, b| b.created_at.cmp(&a.created_at)); // newest first
    existing
        .into_iter()
        .skip(keep_last)
        .filter(|b| now - b.created_at >= min_age)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::Duration;

    #[test]
    fn path_for_strips_traversal_and_empty_segments() {
        let store = FilesystemObjectStore::new("/var/backups");
        assert_eq!(store.path_for("a//b/../c"), PathBuf::from("/var/backups/a/b/c"));
        assert_eq!(store.path_for("catalogue/20260617T031500Z/manifest.json"),
                   PathBuf::from("/var/backups/catalogue/20260617T031500Z/manifest.json"));
    }

    #[test]
    fn filesystem_store_round_trips_and_lists() {
        let root = std::env::temp_dir().join(format!("cerebro-bk-{}", uuid::Uuid::new_v4()));
        let store = FilesystemObjectStore::new(&root);
        store.put("catalogue/a/manifest.json", b"{}").unwrap();
        store.put("catalogue/a/catalogue.archive.gz", b"data").unwrap();
        store.put("other/x", b"z").unwrap();

        assert_eq!(store.get("catalogue/a/catalogue.archive.gz").unwrap(), b"data");
        let mut listed = store.list("catalogue/").unwrap();
        listed.sort();
        assert_eq!(listed, vec![
            "catalogue/a/catalogue.archive.gz".to_string(),
            "catalogue/a/manifest.json".to_string(),
        ]);
        store.delete("catalogue/a/manifest.json").unwrap();
        store.delete("catalogue/a/manifest.json").unwrap(); // idempotent
        assert!(store.list("catalogue/a/manifest.json").unwrap().is_empty());

        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn parse_refs_recovers_ids_and_times() {
        let keys = vec![
            "catalogue/20260617T031500Z/manifest.json".to_string(),
            "catalogue/20260617T031500Z/catalogue.archive.gz".to_string(), // not a manifest
            "catalogue/garbage/manifest.json".to_string(),                 // unparseable id
        ];
        let refs = parse_backup_refs("catalogue", &keys);
        assert_eq!(refs.len(), 1);
        assert_eq!(refs[0].backup_id, "20260617T031500Z");
    }

    #[test]
    fn retention_keeps_last_n_then_honours_age_floor() {
        let now = Utc::now();
        let mk = |days: i64| BackupRef {
            backup_id: format!("b{days}"),
            created_at: now - Duration::days(days),
        };
        // 0,1,2,3,10 days old. keep_last=2 -> keep {0,1}; candidates {2,3,10};
        // min_age=7d -> only the 10-day-old is actually deleted.
        let existing = vec![mk(3), mk(0), mk(10), mk(1), mk(2)];
        let del = select_for_deletion(existing, 2, Duration::days(7), now);
        assert_eq!(del.len(), 1);
        assert_eq!(del[0].backup_id, "b10");
    }

    #[test]
    fn manifest_round_trips_through_json() {
        let m = BackupManifest {
            format_version: MANIFEST_FORMAT_VERSION,
            backup_id: "20260617T031500Z".to_string(),
            created_at: Utc::now(),
            mongo_db: "cerebro".to_string(),
            archive_key: "catalogue/20260617T031500Z/catalogue.archive.gz".to_string(),
            archive_bytes: 4096,
            archive_blake3: "deadbeef".to_string(),
            audit_chain_verified: Some(true),
            audit_event_count: 42,
            cerebro_version: "0.0.0".to_string(),
        };
        let json = serde_json::to_vec(&m).unwrap();
        let back: BackupManifest = serde_json::from_slice(&json).unwrap();
        assert_eq!(m, back);
    }
}

/// S3-compatible object-store backend (S4-4, D1) — AWS S3 / MinIO / SeaweedFS-S3 /
/// Glacier. Feature-gated behind `s3` so the default build pulls in no S3 SDK and
/// the cold tier defaults to [`FilesystemObjectStore`].
///
/// NOTE: this targets the `rust-s3` *blocking* API, whose method names/signatures
/// shift across minor versions. It is the one module that cannot be compile-checked
/// here (the optional dependency is absent by default). When enabling
/// `--features s3`, verify this against the `rust-s3` version pinned in Cargo.toml.
#[cfg(feature = "s3")]
pub struct S3ObjectStore {
    bucket: Box<s3::bucket::Bucket>,
    prefix: String,
}

#[cfg(feature = "s3")]
impl S3ObjectStore {
    /// Open a bucket handle for an S3-compatible endpoint. `prefix` is prepended
    /// to every key so multiple stores can share a bucket.
    pub fn new(
        endpoint: &str,
        region: &str,
        bucket: &str,
        access_key: &str,
        secret_key: &str,
        prefix: &str,
    ) -> anyhow::Result<Self> {
        let region = s3::Region::Custom { region: region.to_string(), endpoint: endpoint.to_string() };
        let creds = s3::creds::Credentials::new(Some(access_key), Some(secret_key), None, None, None)?;
        // Path-style addressing works with MinIO/SeaweedFS-S3 and AWS alike.
        let bucket = s3::bucket::Bucket::new(bucket, region, creds)?.with_path_style();
        Ok(Self { bucket, prefix: prefix.trim_end_matches('/').to_string() })
    }

    fn key(&self, k: &str) -> String {
        if self.prefix.is_empty() {
            k.to_string()
        } else {
            format!("{}/{}", self.prefix, k)
        }
    }
}

#[cfg(feature = "s3")]
impl ObjectStore for S3ObjectStore {
    fn put(&self, key: &str, bytes: &[u8]) -> anyhow::Result<()> {
        self.bucket.put_object_blocking(self.key(key), bytes)?;
        Ok(())
    }

    fn get(&self, key: &str) -> anyhow::Result<Vec<u8>> {
        let resp = self.bucket.get_object_blocking(self.key(key))?;
        Ok(resp.bytes().to_vec())
    }

    fn list(&self, prefix: &str) -> anyhow::Result<Vec<String>> {
        let pages = self.bucket.list_blocking(self.key(prefix), None)?;
        let mut keys = Vec::new();
        for (page, _status) in pages {
            for obj in page.contents {
                keys.push(obj.key);
            }
        }
        Ok(keys)
    }

    fn delete(&self, key: &str) -> anyhow::Result<()> {
        self.bucket.delete_object_blocking(self.key(key))?;
        Ok(())
    }
}
