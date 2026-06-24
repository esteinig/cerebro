//! CIQA persistence records: QC datasets and versioned baselines.
//!
//! Generalises the training dataset pattern (`TrainingPrefetchRecord` — `PrefetchData` in GridFS)
//! into a QC dataset that additionally carries **reference truth** for sensitivity/specificity, and
//! adds an **immutable, versioned baseline** a run is regressed against. New fields are
//! `#[serde(default)]`; baselines are append-only and never overwritten (D12).

use serde::{Deserialize, Serialize};

use crate::api::cerebro::schema::{SampleType, TestResult};
use crate::api::ciqa::schema::MetaGptRunManifest;
use crate::api::diagnostics::eval::DiagnosticStatistics;
use crate::api::files::model::SeaweedFile;

/// Pass thresholds for a QC dataset/baseline. Values are **fractions** (`0.0..=1.0`), matching
/// `DiagnosticStatistics`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CiqaThresholds {
    pub sensitivity: f64,
    pub specificity: f64,
}

/// Staged QC dataset descriptor (FASTQ staged under a `ciqa-<run_id>` prefix, as today).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CiqaDataset {
    pub files: Vec<SeaweedFile>,
    pub description: String,
    pub version: String,
    pub qa_thresholds: CiqaThresholds,
    pub dataset_prefix: String,
}

/// One QC dataset record: a `PrefetchData` in GridFS (by `id`, as `TrainingPrefetchRecord`) plus
/// the reference truth needed to score sensitivity/specificity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityControlDatasetRecord {
    /// GridFS id of the stored `PrefetchData` (mirrors `TrainingPrefetchRecord`).
    pub id: String,
    /// Dataset collection/name.
    pub dataset: String,
    pub description: String,
    /// Dataset version (immutable once a baseline references it).
    pub version: String,
    #[serde(default)]
    pub sample_name: Option<String>,
    #[serde(default)]
    pub sample_type: Option<SampleType>,
    /// Reference result truth for this sample.
    #[serde(default)]
    pub reference_result: Option<TestResult>,
    /// Normalised `s__...` reference candidate truth.
    #[serde(default)]
    pub reference_candidates: Option<Vec<String>>,
    #[serde(default)]
    pub exclude_lod: Option<bool>,
}

/// An immutable, versioned baseline: the frozen statistics a run is regressed against, plus the
/// manifest that produced it and the thresholds to pass. Append-only; promotion is operator-gated.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityControlBaseline {
    pub id: String,
    pub dataset: String,
    pub dataset_version: String,
    /// What produced the baseline.
    pub model_manifest: MetaGptRunManifest,
    /// Frozen sensitivity/specificity/PPV/NPV + per-sample decisions.
    pub statistics: DiagnosticStatistics,
    /// Floor a run must meet to pass.
    pub thresholds: CiqaThresholds,
    pub created: String,
    pub created_by: String,
    /// Integrity hash (SKILLS.md §4).
    pub content_hash: String,
    /// Operator-gated promotion to the active baseline for the dataset.
    #[serde(default)]
    pub promoted: bool,
}

impl QualityControlBaseline {
    pub fn to_json<P: AsRef<std::path::Path>>(
        &self,
        path: P,
    ) -> Result<(), crate::api::cerebro::model::ModelError> {
        let data = serde_json::to_string_pretty(self)
            .map_err(crate::api::cerebro::model::ModelError::JsonSerialization)?;
        std::fs::write(path, data)?;
        Ok(())
    }
    pub fn from_json<P: AsRef<std::path::Path>>(
        path: P,
    ) -> Result<Self, crate::api::cerebro::model::ModelError> {
        let data = std::fs::read_to_string(path)?;
        serde_json::from_str(&data)
            .map_err(crate::api::cerebro::model::ModelError::JsonDeserialization)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dataset_record_round_trips_with_defaults() {
        let json = r#"{"id":"gridfs1","dataset":"qc-v1","description":"core panel","version":"1.0.0",
            "reference_result":"positive","reference_candidates":["s__Escherichia coli"]}"#;
        let r: QualityControlDatasetRecord = serde_json::from_str(json).unwrap();
        assert_eq!(r.id, "gridfs1");
        assert_eq!(r.reference_result, Some(TestResult::Positive));
        assert_eq!(r.sample_type, None); // defaulted
        assert_eq!(r.exclude_lod, None);
    }
}
