//! Run provenance manifest.
//!
//! A [`RunManifest`] is the reproducibility / chain-of-custody record for a
//! pipeline run: which pipeline (name + version) produced which outputs, under
//! what parameters, with which tool and reference-database versions, and the
//! BLAKE3 hashes of the input and output artefacts. It is itself captured into
//! cerebro-fs as a [`FileType::RunManifest`](super::model::FileType) artefact.
//!
//! The manifest is **sealed** with a BLAKE3 content hash over its canonical body
//! (every field except the seal itself), making it self-verifying / tamper
//! evident. Computing and verifying that seal lives in the `cerebro-fs` crate
//! (where BLAKE3 is), so this module stays a pure, dependency-light data model.
//!
//! Note: the seal provides integrity, not non-repudiation. Cryptographic signing
//! (e.g. ed25519) would be an additive follow-on requiring a signing dependency.

use std::collections::BTreeMap;

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

use super::model::FileType;

/// A tool and the version used in the run.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct ToolVersion {
    pub name: String,
    pub version: String,
}

/// A reference database and the version (and optional hash) used in the run.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct ReferenceDb {
    pub name: String,
    pub version: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hash: Option<String>,
}

/// An input or output artefact referenced by the manifest, by BLAKE3 hash.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ManifestArtefact {
    pub name: String,
    pub hash: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ftype: Option<FileType>,
}

/// Provenance metadata supplied by the pipeline run.
///
/// Deserialisable from a JSON file so a pipeline can emit it directly; the
/// remaining manifest fields (inputs/outputs/identifiers/seal) are filled in when
/// the manifest is built from the registered artefacts.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ManifestProvenance {
    pub pipeline_name: String,
    pub pipeline_version: String,
    #[serde(default)]
    pub parameters: BTreeMap<String, String>,
    #[serde(default)]
    pub tool_versions: Vec<ToolVersion>,
    #[serde(default)]
    pub reference_dbs: Vec<ReferenceDb>,
}

/// The full provenance manifest for a run/sample.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunManifest {
    pub run_id: Option<String>,
    pub sample_id: Option<String>,
    pub pipeline_id: Option<String>,
    pub pipeline_name: String,
    pub pipeline_version: String,
    pub created_at: DateTime<Utc>,
    /// Pipeline parameters (sorted for a deterministic seal).
    pub parameters: BTreeMap<String, String>,
    pub tool_versions: Vec<ToolVersion>,
    pub reference_dbs: Vec<ReferenceDb>,
    pub inputs: Vec<ManifestArtefact>,
    pub outputs: Vec<ManifestArtefact>,
    /// BLAKE3 seal over the canonical body (all other fields). `None` until
    /// sealed. Computed/verified by the `cerebro-fs` crate.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub content_hash: Option<String>,
}

impl RunManifest {
    /// Create an unsealed manifest from provenance metadata and identifiers;
    /// `inputs`/`outputs` start empty and are populated when built from the
    /// registered artefacts.
    pub fn new(
        run_id: Option<String>,
        sample_id: Option<String>,
        pipeline_id: Option<String>,
        provenance: ManifestProvenance,
        created_at: DateTime<Utc>,
    ) -> Self {
        Self {
            run_id,
            sample_id,
            pipeline_id,
            pipeline_name: provenance.pipeline_name,
            pipeline_version: provenance.pipeline_version,
            created_at,
            parameters: provenance.parameters,
            tool_versions: provenance.tool_versions,
            reference_dbs: provenance.reference_dbs,
            inputs: Vec::new(),
            outputs: Vec::new(),
            content_hash: None,
        }
    }

    /// A clone with the seal cleared — the canonical object the content hash is
    /// computed over (so the hash never covers itself).
    pub fn body_for_seal(&self) -> Self {
        let mut body = self.clone();
        body.content_hash = None;
        body
    }

    /// Record the computed seal.
    pub fn set_content_hash(&mut self, hash: String) {
        self.content_hash = Some(hash);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn body_for_seal_excludes_the_seal() {
        let mut m = RunManifest::new(
            Some("RUN".into()),
            None,
            None,
            ManifestProvenance {
                pipeline_name: "cerebro".into(),
                pipeline_version: "1.0".into(),
                ..Default::default()
            },
            Utc::now(),
        );
        m.set_content_hash("deadbeef".into());
        assert!(m.content_hash.is_some());
        assert!(m.body_for_seal().content_hash.is_none());
        // Setting the seal must not perturb the rest of the body.
        let mut cleared = m.clone();
        cleared.content_hash = None;
        assert_eq!(
            serde_json::to_value(&cleared).unwrap(),
            serde_json::to_value(&m.body_for_seal()).unwrap()
        );
    }
}
