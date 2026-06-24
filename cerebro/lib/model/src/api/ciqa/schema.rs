//! CIQA run attribution: the typed manifest the diagnose step emits and the regression step
//! consumes.
//!
//! This replaces directory-name regex (`parse_dir_components`) on the production path: the model
//! id, quantization, parameters, clinical flag, and replicate are now **carried explicitly and
//! correctly** rather than recovered by pattern-matching folder names. The manifest also records,
//! per sample, the prefetch/result paths and an integrity hash of the output (SKILLS.md §4).

use serde::{Deserialize, Serialize};

use crate::api::cerebro::model::ModelError;

/// Integrity hash of arbitrary content (blake3 hex), used for manifest/result hashes
/// (SKILLS.md §4). Consistent with the model crate's existing blake3-based hashing.
pub fn content_hash(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

// ── CIQA API request / query schemas (Stage 12). Mirrors training::schema. ───────────────

use crate::api::cerebro::schema::{PrefetchData, SampleType, TestResult};
use crate::api::diagnostics::eval::DiagnosticStatistics;

/// POST /ciqa/dataset — store one QC dataset sample: the `PrefetchData` goes to GridFS (by `id`),
/// the reference truth + metadata become the record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateCiqaDataset {
    pub id: String,
    pub dataset: String,
    pub description: String,
    pub version: String,
    pub prefetch: PrefetchData,
    #[serde(default)]
    pub sample_name: Option<String>,
    #[serde(default)]
    pub sample_type: Option<SampleType>,
    #[serde(default)]
    pub reference_result: Option<TestResult>,
    #[serde(default)]
    pub reference_candidates: Option<Vec<String>>,
    #[serde(default)]
    pub exclude_lod: Option<bool>,
}

/// POST /ciqa/baseline — register an immutable baseline. The server sets `created`,
/// `content_hash`, and `promoted=false`; promotion is a separate gated step.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateCiqaBaseline {
    pub id: String,
    pub dataset: String,
    pub dataset_version: String,
    pub model_manifest: MetaGptRunManifest,
    pub statistics: DiagnosticStatistics,
    pub thresholds: super::model::CiqaThresholds,
    pub created_by: String,
}

/// POST /ciqa/regression — persist a run's regression outcome against a baseline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateCiqaRegression {
    pub run_id: String,
    pub baseline_id: String,
    pub manifest: MetaGptRunManifest,
    pub statistics: DiagnosticStatistics,
    pub regressed: bool,
    #[serde(default)]
    pub created: Option<String>,
}

/// Query filter for dataset/baseline listings (team-scoped server-side).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryCiqa {
    #[serde(default)]
    pub dataset: Option<String>,
    #[serde(default)]
    pub version: Option<String>,
}


/// What produced a META-GPT run: model/config attribution + per-sample outputs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetaGptRunManifest {
    pub run_id: String,
    pub model_id: String,
    #[serde(default)]
    pub quantization: Option<String>,
    #[serde(default)]
    pub params: Option<String>,
    pub clinical: bool,
    #[serde(default)]
    pub replicate: Option<u32>,
    /// Hash of the `MetaGpConfig`/`TieredFilterConfig` used to produce the prefetch.
    pub config_hash: String,
    #[serde(default)]
    pub samples: Vec<MetaGptRunSample>,
    #[serde(default)]
    pub created: Option<String>,
}

/// One sample's outputs within a run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetaGptRunSample {
    pub sample_id: String,
    pub prefetch_path: String,
    /// `{sample}.model.json` — the META-GPT `DiagnosticResult` for this sample.
    pub result_path: String,
    /// Integrity hash of the result file (SKILLS.md §4).
    pub content_hash: String,
}

impl MetaGptRunManifest {
    pub fn to_json<P: AsRef<std::path::Path>>(&self, path: P) -> Result<(), ModelError> {
        let data =
            serde_json::to_string_pretty(self).map_err(ModelError::JsonSerialization)?;
        std::fs::write(path, data)?;
        Ok(())
    }
    pub fn from_json<P: AsRef<std::path::Path>>(path: P) -> Result<Self, ModelError> {
        let data = std::fs::read_to_string(path)?;
        serde_json::from_str(&data).map_err(ModelError::JsonDeserialization)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn content_hash_is_deterministic_and_distinct() {
        let a = content_hash(b"hello");
        assert_eq!(a, content_hash(b"hello"));
        assert_ne!(a, content_hash(b"world"));
        assert_eq!(a.len(), 64); // blake3 256-bit hex
    }

    #[test]
    fn manifest_round_trips_with_defaults() {
        let json = r#"{
            "run_id":"RUN1","model_id":"qwen3-14b","clinical":true,"config_hash":"abc",
            "samples":[{"sample_id":"S1","prefetch_path":"S1.prefetch.json","result_path":"S1.model.json","content_hash":"h1"}]
        }"#;
        let m: MetaGptRunManifest = serde_json::from_str(json).unwrap();
        assert_eq!(m.run_id, "RUN1");
        assert_eq!(m.quantization, None); // defaulted
        assert_eq!(m.samples.len(), 1);
        let again: MetaGptRunManifest =
            serde_json::from_str(&serde_json::to_string(&m).unwrap()).unwrap();
        assert_eq!(again.model_id, "qwen3-14b");
        assert_eq!(again.clinical, true);
    }
}
