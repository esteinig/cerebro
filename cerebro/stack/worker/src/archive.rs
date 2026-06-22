//! Archival engine (S4-4): copy a file's bytes to the cold object store so the
//! cold tier is *real* rather than just a catalogue label.
//!
//! Backend-agnostic — it works against any [`ObjectStore`]: a mounted cold disk
//! or NFS via [`crate::backup::FilesystemObjectStore`] (the compile-safe default),
//! or S3 / Glacier via the feature-gated S3 backend. After the bytes are copied,
//! the caller repoints the catalogue through the dedicated relocate endpoint (D3),
//! which flips `archived = true` and records the `archive_key`.

use cerebro_fs::client::FileSystemClient;

use crate::backup::ObjectStore;

/// Outcome of archiving one object.
#[derive(Debug, Clone)]
pub struct ArchiveOutcome {
    pub archive_key: String,
    pub bytes: usize,
}

/// Derive the stable cold-store key for a file's bytes: `<prefix>/<team>/<file_id>`.
///
/// Namespaced by team and id so keys are unique and human-traceable, and
/// sanitised so they are always valid object keys. Deterministic, so re-archiving
/// the same file targets the same key (idempotent on the store side). Pure.
pub fn archive_key_for(prefix: &str, team: &str, file_id: &str) -> String {
    let p = prefix.trim_matches('/');
    let t = sanitize(team);
    let f = sanitize(file_id);
    if p.is_empty() {
        format!("{t}/{f}")
    } else {
        format!("{p}/{t}/{f}")
    }
}

fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.') {
                c
            } else {
                '_'
            }
        })
        .collect()
}

/// Copy `fid`'s bytes from SeaweedFS to `store` at `archive_key`.
///
/// Blocking (reqwest under the hood) — call from a blocking context. The caller
/// then repoints the catalogue via the relocate endpoint; only after that repoint
/// is the object officially archived, so a failure here leaves the catalogue
/// untouched and the file directly retrievable.
pub fn archive_object(
    fs: &FileSystemClient,
    store: &dyn ObjectStore,
    fid: &str,
    archive_key: &str,
) -> anyhow::Result<ArchiveOutcome> {
    let bytes = fs
        .read_object(fid)
        .map_err(|e| anyhow::anyhow!("read object {fid} from store: {e}"))?;
    store.put(archive_key, &bytes)?;
    Ok(ArchiveOutcome {
        archive_key: archive_key.to_string(),
        bytes: bytes.len(),
    })
}

/// Outcome of restoring one object from the cold store (S4-5).
#[derive(Debug, Clone)]
pub struct RestoreOutcome {
    /// `Some(fid)` when a fresh weed object was written (fid-addressed file — the
    /// caller must repoint `fid`); `None` when a path-addressed filer object was
    /// overwritten **in place** (the path is unchanged, so no fid update is needed).
    pub new_fid: Option<String>,
    pub bytes: usize,
    pub hash_verified: bool,
}

/// Restore an archived object from `store` back into SeaweedFS (S4-5; in-place for
/// path-addressed files in S4-6 H4).
///
/// Fetches the bytes at `archive_key`, verifies them against `expected_hash` (the
/// catalogue BLAKE3) when supplied — so a corrupt cold copy can never silently
/// replace a good catalogue entry — then writes them back to the file's
/// **effective location** given in `effective_id`:
///
/// * a **filer path** (contains `/`) is overwritten in place, keeping the path
///   valid and pointing at fresh data — `new_fid` is `None`;
/// * otherwise a **new weed object** is written and its fid returned in `new_fid`.
///
/// The caller repoints the catalogue (relocate): it sets `fid` only when `new_fid`
/// is `Some`, and lands the file on a retrievable tier. Blocking — call from a
/// blocking context.
pub fn restore_object(
    fs: &FileSystemClient,
    store: &dyn ObjectStore,
    archive_key: &str,
    effective_id: &str,
    name: &str,
    expected_hash: Option<&str>,
) -> anyhow::Result<RestoreOutcome> {
    let bytes = store.get(archive_key)?;
    let hash_verified = match expected_hash {
        Some(expected) => {
            let actual = cerebro_fs::hash::hash_bytes(&bytes);
            if actual != expected {
                anyhow::bail!(
                    "restored bytes hash mismatch for {archive_key}: expected {expected}, got {actual}"
                );
            }
            true
        }
        None => false,
    };
    let new_fid = if cerebro_fs::client::is_filer_path(effective_id) {
        // Path-addressed: overwrite the filer object in place; the path stays valid.
        fs.write_object_at_path(effective_id, &bytes)
            .map_err(|e| anyhow::anyhow!("re-materialise {archive_key} at {effective_id}: {e}"))?;
        None
    } else {
        // Fid-addressed: write a fresh weed object; the caller repoints fid.
        let fid = fs
            .write_object(name, &bytes)
            .map_err(|e| anyhow::anyhow!("re-materialise {archive_key}: {e}"))?;
        Some(fid)
    };
    Ok(RestoreOutcome {
        new_fid,
        bytes: bytes.len(),
        hash_verified,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn archive_key_is_namespaced_and_sanitised() {
        assert_eq!(
            archive_key_for("archive", "team-a", "file_1"),
            "archive/team-a/file_1"
        );
        // unsafe characters in team/id are replaced
        assert_eq!(
            archive_key_for("archive", "te/am", "f i d"),
            "archive/te_am/f_i_d"
        );
        // empty prefix drops the leading segment
        assert_eq!(archive_key_for("", "t", "f"), "t/f");
        // surrounding slashes on the prefix are trimmed
        assert_eq!(archive_key_for("/archive/", "t", "f"), "archive/t/f");
    }

    #[test]
    fn archive_key_is_deterministic() {
        let a = archive_key_for("archive", "team", "abc");
        let b = archive_key_for("archive", "team", "abc");
        assert_eq!(a, b);
    }
}
