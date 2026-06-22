//! Run manifest sealing, assembly and capture.
//!
//! Builds a [`RunManifest`] from the artefacts already registered for a run,
//! **seals** it with a BLAKE3 content hash over its canonical body (self-verifying
//! / tamper evident), and captures the sealed manifest into cerebro-fs as a
//! [`FileType::RunManifest`] artefact — closing the provenance loop opened by the
//! output capture in [`crate::capture`].
//!
//! Two hashes coexist and serve different purposes: the manifest's `content_hash`
//! seals the provenance body, while the registered `SeaweedFile.hash` is the
//! storage-integrity hash of the manifest JSON file (what verify/repair use).

use std::path::PathBuf;

use chrono::Utc;
use uuid::Uuid;

use cerebro_model::api::files::manifest::{ManifestArtefact, ManifestProvenance, RunManifest};
use cerebro_model::api::files::model::{FileType, SeaweedFile};
use cerebro_model::api::files::retention::RetentionClass;
use cerebro_model::api::files::schema::RegisterFileSchema;

use crate::client::{FileSystemClient, UploadConfig};
use crate::error::FileSystemError;
use crate::hash::{fast_file_hash, hash_bytes};

/// Serialise the canonical (seal-excluded) manifest body to JSON bytes.
pub fn canonical_bytes(manifest: &RunManifest) -> Result<Vec<u8>, FileSystemError> {
    serde_json::to_vec(&manifest.body_for_seal())
        .map_err(|e| FileSystemError::ManifestSerialization(e.to_string()))
}

/// Seal a manifest in place by computing its BLAKE3 content hash.
pub fn seal(manifest: &mut RunManifest) -> Result<(), FileSystemError> {
    let hash = hash_bytes(&canonical_bytes(manifest)?);
    manifest.set_content_hash(hash);
    Ok(())
}

/// Verify a manifest's seal by recomputing the content hash over its canonical
/// body. Returns `false` when unsealed or when the body has been altered.
pub fn verify(manifest: &RunManifest) -> Result<bool, FileSystemError> {
    match &manifest.content_hash {
        None => Ok(false),
        Some(recorded) => Ok(hash_bytes(&canonical_bytes(manifest)?) == *recorded),
    }
}

/// True for artefacts that are run inputs (raw reads) rather than outputs.
fn is_input(file: &SeaweedFile) -> bool {
    matches!(
        file.ftype,
        Some(FileType::ReadPaired) | Some(FileType::ReadSingle)
    )
}

impl FileSystemClient {
    /// Build a sealed [`RunManifest`] for a run/sample from its registered
    /// artefacts.
    ///
    /// Resolves the run's files via the API, partitions them into inputs (raw
    /// reads) and outputs (everything else, excluding any existing manifest), and
    /// records each as a [`ManifestArtefact`] referencing its stored BLAKE3 hash
    /// and id. The manifest is sealed before being returned.
    pub fn build_run_manifest(
        &self,
        run_id: Option<String>,
        sample_id: Option<String>,
        pipeline_id: Option<String>,
        provenance: ManifestProvenance,
    ) -> Result<RunManifest, FileSystemError> {
        let files = self
            .api_client
            .list_files(run_id.clone(), None, 0, 100_000, false)?;
        let files: Vec<SeaweedFile> = match &sample_id {
            Some(sid) => files
                .into_iter()
                .filter(|f| f.sample_id.as_deref() == Some(sid.as_str()))
                .collect(),
            None => files,
        };

        let to_artefact = |f: &SeaweedFile| ManifestArtefact {
            name: f.name.clone(),
            hash: f.hash.clone(),
            file_id: Some(f.id.clone()),
            ftype: f.ftype.clone(),
        };

        let mut manifest = RunManifest::new(run_id, sample_id, pipeline_id, provenance, Utc::now());
        for file in &files {
            // Never reference a previous manifest as an output of this run.
            if matches!(file.ftype, Some(FileType::RunManifest)) {
                continue;
            }
            if is_input(file) {
                manifest.inputs.push(to_artefact(file));
            } else {
                manifest.outputs.push(to_artefact(file));
            }
        }

        seal(&mut manifest)?;
        Ok(manifest)
    }

    /// Capture a sealed manifest into cerebro-fs as a `RunManifest` artefact.
    ///
    /// Writes the manifest JSON to a temporary file, uploads it via the configured
    /// access mode, and registers it (diagnostic retention, hot tier) linked to
    /// the run/sample. Returns the registered file id.
    pub fn capture_manifest(
        &self,
        manifest: &RunManifest,
        upload_config: &UploadConfig,
    ) -> Result<String, FileSystemError> {
        let json = serde_json::to_vec_pretty(manifest)
            .map_err(|e| FileSystemError::ManifestSerialization(e.to_string()))?;

        let name = match (&manifest.run_id, &manifest.sample_id) {
            (Some(run), Some(sample)) => format!("{run}.{sample}.manifest.json"),
            (Some(run), None) => format!("{run}.manifest.json"),
            _ => "run.manifest.json".to_string(),
        };

        let tmp = std::env::temp_dir().join(format!("cerebro-manifest-{}", std::process::id()));
        std::fs::create_dir_all(&tmp)?;
        let path: PathBuf = tmp.join(&name);
        std::fs::write(&path, &json)?;

        let file_hash = fast_file_hash(&path)?;
        let stored = self.store_object(
            &path,
            manifest.run_id.as_deref(),
            manifest.sample_id.as_deref(),
            upload_config,
        )?;

        let now = Utc::now();
        let id = Uuid::new_v4().to_string();
        let file_schema = RegisterFileSchema {
            id: id.clone(),
            run_id: manifest.run_id.clone(),
            sample_id: manifest.sample_id.clone(),
            pipeline_id: manifest.pipeline_id.clone(),
            description: Some("run provenance manifest".to_string()),
            date: now.to_string(),
            name: stored.name,
            hash: file_hash,
            fid: stored.fid,
            ftype: Some(FileType::RunManifest),
            size: stored.size,
            watcher: None,
            path: stored.path,
            tier: upload_config.tier,
            retention: RetentionClass::Diagnostic,
            retain_until: upload_config
                .retention_policy
                .retain_until(RetentionClass::Diagnostic, now),
            legal_hold: upload_config.legal_hold,
            replicas: None,
            archived: false,
            reported_at: None,
        };

        self.api_client.register_file(file_schema)?;

        // Best-effort cleanup of the temporary manifest file.
        let _ = std::fs::remove_dir_all(&tmp);

        Ok(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cerebro_model::api::files::manifest::ManifestProvenance;

    fn sample_manifest() -> RunManifest {
        let mut m = RunManifest::new(
            Some("RUN01".into()),
            Some("S1".into()),
            None,
            ManifestProvenance {
                pipeline_name: "cerebro".into(),
                pipeline_version: "1.0.0".into(),
                ..Default::default()
            },
            Utc::now(),
        );
        m.outputs.push(ManifestArtefact {
            name: "S1.cerebro.json".into(),
            hash: "abc".into(),
            file_id: Some("id1".into()),
            ftype: Some(FileType::CerebroModel),
        });
        m
    }

    #[test]
    fn seal_then_verify_roundtrips() {
        let mut m = sample_manifest();
        assert!(!verify(&m).unwrap()); // unsealed
        seal(&mut m).unwrap();
        assert!(m.content_hash.is_some());
        assert!(verify(&m).unwrap());
    }

    #[test]
    fn tampering_breaks_the_seal() {
        let mut m = sample_manifest();
        seal(&mut m).unwrap();
        // Alter the body after sealing.
        m.outputs[0].hash = "tampered".into();
        assert!(!verify(&m).unwrap());
    }

    #[test]
    fn seal_is_stable_and_excludes_itself() {
        let mut m = sample_manifest();
        seal(&mut m).unwrap();
        let first = m.content_hash.clone();
        // Re-sealing the same body yields the same hash (seal excludes itself).
        seal(&mut m).unwrap();
        assert_eq!(first, m.content_hash);
    }
}
