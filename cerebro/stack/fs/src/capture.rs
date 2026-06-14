//! Pipeline output capture for Cerebro FS (S2-2).
//!
//! After a pipeline run completes, its diagnostic artefacts must be captured back
//! into cerebro-fs — hashed, typed, and stamped with a tier and retention class —
//! and linked to their run/sample. This module walks a run's output directory,
//! **classifies** each file into a [`FileType`] and [`RetentionClass`] via a
//! tunable ruleset, then uploads and registers it on the existing store path.
//!
//! Capture is sequential within a node and continues past individual failures
//! (recording them), so one bad file never loses the rest of a run. Running
//! capture **execution-node-parallel** across the cluster is the orchestration
//! layer's job (S2-5), not in-process concurrency here.
//!
//! The default ruleset ([`CaptureRule::default_ruleset`]) is a sensible starting
//! point keyed off common output conventions; **tune it to your pipeline's actual
//! output layout** — it is matched by lowercased substring, first rule wins.

use std::path::{Path, PathBuf};

use chrono::Utc;
use uuid::Uuid;
use serde::Deserialize;

use cerebro_model::api::files::model::FileType;
use cerebro_model::api::files::retention::RetentionClass;
use cerebro_model::api::files::schema::RegisterFileSchema;

use crate::client::{FileSystemClient, UploadConfig};
use crate::error::FileSystemError;
use crate::hash::fast_file_hash;

/// How a matched file should be handled.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ArtefactClass {
    /// Capture the file under the given type and retention class.
    Capture { ftype: Option<FileType>, retention: RetentionClass },
    /// Skip the file entirely (e.g. scratch, logs, Nextflow work dirs).
    Ignore,
}

/// A single classification rule: if `pattern` is a (lowercased) substring of a
/// file's path relative to the output root, it is handled as `class`.
#[derive(Debug, Clone)]
pub struct CaptureRule {
    pub pattern: String,
    pub class: ArtefactClass,
}
impl CaptureRule {
    fn capture(pattern: &str, ftype: FileType, retention: RetentionClass) -> Self {
        Self { pattern: pattern.to_string(), class: ArtefactClass::Capture { ftype: Some(ftype), retention } }
    }
    fn ignore(pattern: &str) -> Self {
        Self { pattern: pattern.to_string(), class: ArtefactClass::Ignore }
    }

    /// A default ruleset for the Cerebro pipeline outputs. Ordered — first match
    /// wins. **Tune to your actual output layout.**
    ///
    /// Ignores Nextflow scratch/logs; types QC as intermediate; pathogen,
    /// panviral, consensus and the Cerebro model as diagnostic.
    pub fn default_ruleset() -> Vec<CaptureRule> {
        vec![
            // Scratch / logs — never captured.
            CaptureRule::ignore("work/"),
            CaptureRule::ignore(".command."),
            CaptureRule::ignore(".nextflow"),
            CaptureRule::ignore(".log"),
            // Primary result document.
            CaptureRule::capture(".cerebro.json", FileType::CerebroModel, RetentionClass::Diagnostic),
            // Diagnostic outputs.
            CaptureRule::capture("consensus", FileType::Consensus, RetentionClass::Diagnostic),
            CaptureRule::capture("panviral", FileType::PanviralOutput, RetentionClass::Diagnostic),
            CaptureRule::capture("pathogen", FileType::PathogenOutput, RetentionClass::Diagnostic),
            // Quality control — useful for re-inspection but intermediate.
            CaptureRule::capture("quality", FileType::QualityOutput, RetentionClass::Intermediate),
            CaptureRule::capture("fastp", FileType::QualityOutput, RetentionClass::Intermediate),
            CaptureRule::capture("fastqc", FileType::QualityOutput, RetentionClass::Intermediate),
        ]
    }
}

/// Classify a relative path against `rules`, falling back to
/// `Other`/`Intermediate` when no rule matches.
pub fn classify(relative_path: &str, rules: &[CaptureRule]) -> ArtefactClass {
    let haystack = relative_path.to_lowercase();
    for rule in rules {
        if haystack.contains(&rule.pattern) {
            return rule.class.clone();
        }
    }
    ArtefactClass::Capture { ftype: Some(FileType::Other), retention: RetentionClass::Intermediate }
}

/// A single rule as declared by the pipeline in `outputs.json` (intermission-3).
#[derive(Debug, Deserialize)]
struct PipelineRule {
    pattern: String,
    #[serde(default)]
    ignore: bool,
    #[serde(default)]
    ftype: Option<FileType>,
    #[serde(default)]
    retention: RetentionClass,
}

/// The pipeline-emitted outputs manifest: an authoritative classification ruleset
/// for the run's artefacts. Unknown fields (e.g. `report_out`) are ignored.
#[derive(Debug, Deserialize)]
struct OutputsManifest {
    #[serde(default)]
    rules: Vec<PipelineRule>,
}

impl CaptureRule {
    /// Load pipeline-declared capture rules from an `outputs.json` (intermission-3).
    ///
    /// The pipeline knows what each artefact *is*, so its declared rules take
    /// precedence over the filename heuristics in [`default_ruleset`]. Returns an
    /// empty vec when the file is absent or invalid (capture then falls back to the
    /// defaults). Patterns are lowercased to match [`classify`].
    pub fn from_manifest(path: &Path) -> Vec<CaptureRule> {
        let Ok(bytes) = std::fs::read(path) else { return Vec::new(); };
        let Ok(manifest) = serde_json::from_slice::<OutputsManifest>(&bytes) else { return Vec::new(); };
        manifest
            .rules
            .into_iter()
            .map(|r| {
                let class = if r.ignore {
                    ArtefactClass::Ignore
                } else {
                    ArtefactClass::Capture { ftype: r.ftype, retention: r.retention }
                };
                CaptureRule { pattern: r.pattern.to_lowercase(), class }
            })
            .collect()
    }
}

/// Locate a pipeline `outputs.json` for a capture target: the dir itself or its
/// `output/` child (where the pipeline publishes it).
fn find_outputs_manifest(dir: &Path) -> Option<PathBuf> {
    for candidate in [dir.join("outputs.json"), dir.join("output").join("outputs.json")] {
        if candidate.is_file() {
            return Some(candidate);
        }
    }
    None
}

/// Outcome of handling a single file during capture.
#[derive(Debug, Clone)]
pub enum CaptureStatus {
    Captured { ftype: Option<FileType>, retention: RetentionClass },
    Ignored,
    Failed(String),
}

/// Per-file capture outcome.
#[derive(Debug, Clone)]
pub struct CaptureOutcome {
    /// Path relative to the output root.
    pub relative_path: String,
    pub status: CaptureStatus,
}

/// Aggregate report of a capture sweep.
#[derive(Debug, Clone)]
pub struct CaptureReport {
    pub outcomes: Vec<CaptureOutcome>,
}
impl CaptureReport {
    pub fn captured(&self) -> usize {
        self.outcomes.iter().filter(|o| matches!(o.status, CaptureStatus::Captured { .. })).count()
    }
    pub fn ignored(&self) -> usize {
        self.outcomes.iter().filter(|o| matches!(o.status, CaptureStatus::Ignored)).count()
    }
    pub fn failed(&self) -> usize {
        self.outcomes.iter().filter(|o| matches!(o.status, CaptureStatus::Failed(_))).count()
    }
    /// True when nothing failed.
    pub fn ok(&self) -> bool {
        self.failed() == 0
    }
}

/// Recursively collect regular files under `dir`.
fn collect_files(dir: &Path, out: &mut Vec<PathBuf>) -> Result<(), FileSystemError> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        let file_type = entry.file_type()?;
        if file_type.is_dir() {
            collect_files(&path, out)?;
        } else if file_type.is_file() {
            out.push(path);
        }
        // Symlinks are intentionally not followed.
    }
    Ok(())
}

impl FileSystemClient {
    /// Capture a completed run's output directory into cerebro-fs.
    ///
    /// Walks `output_dir`, classifies each file with `rules`
    /// ([`CaptureRule::default_ruleset`] is a good default), and for every
    /// captured file: hashes it (BLAKE3), uploads it via the configured access
    /// mode, and registers a `SeaweedFile` stamped with the classified
    /// [`FileType`]/[`RetentionClass`], the run/sample/pipeline identifiers, the
    /// hot tier, and a `retain_until` computed from `upload_config`'s policy.
    ///
    /// Continues past per-file failures, recording them in the returned
    /// [`CaptureReport`]; inspect [`CaptureReport::ok`].
    pub fn capture_outputs(
        &self,
        run_id: Option<String>,
        sample_id: Option<String>,
        pipeline_id: Option<String>,
        output_dir: &Path,
        rules: &[CaptureRule],
        upload_config: &UploadConfig,
    ) -> Result<CaptureReport, FileSystemError> {

        if !output_dir.is_dir() {
            return Err(FileSystemError::FileDoesNotExist(output_dir.to_path_buf()));
        }

        // Prefer the pipeline's declared classification (intermission-3): if the
        // run output carries an outputs.json, its rules take precedence over the
        // caller's defaults (first match wins in `classify`).
        let effective_rules: Vec<CaptureRule> = match find_outputs_manifest(output_dir) {
            Some(manifest) => {
                let pipeline_rules = CaptureRule::from_manifest(&manifest);
                log::info!("Using {} pipeline-declared capture rule(s) from {}", pipeline_rules.len(), manifest.display());
                pipeline_rules.into_iter().chain(rules.iter().cloned()).collect()
            }
            None => rules.to_vec(),
        };

        let mut files = Vec::new();
        collect_files(output_dir, &mut files)?;

        let mut outcomes = Vec::new();
        for file in &files {
            let relative_path = file
                .strip_prefix(output_dir)
                .unwrap_or(file)
                .to_string_lossy()
                .to_string();

            let (ftype, retention) = match classify(&relative_path, &effective_rules) {
                ArtefactClass::Ignore => {
                    outcomes.push(CaptureOutcome { relative_path, status: CaptureStatus::Ignored });
                    continue;
                }
                ArtefactClass::Capture { ftype, retention } => (ftype, retention),
            };

            match self.capture_one(file, run_id.as_deref(), sample_id.as_deref(), pipeline_id.as_deref(), ftype.clone(), retention, upload_config) {
                Ok(()) => outcomes.push(CaptureOutcome { relative_path, status: CaptureStatus::Captured { ftype, retention } }),
                Err(err) => outcomes.push(CaptureOutcome { relative_path, status: CaptureStatus::Failed(err.to_string()) }),
            }
        }

        Ok(CaptureReport { outcomes })
    }

    /// Hash, upload and register a single captured output file.
    fn capture_one(
        &self,
        file: &PathBuf,
        run_id: Option<&str>,
        sample_id: Option<&str>,
        pipeline_id: Option<&str>,
        ftype: Option<FileType>,
        retention: RetentionClass,
        upload_config: &UploadConfig,
    ) -> Result<(), FileSystemError> {
        let file_hash = fast_file_hash(file)?;
        let stored = self.store_object(file, run_id, sample_id, upload_config)?;

        let now = Utc::now();
        let retain_until = upload_config.retention_policy.retain_until(retention, now);

        let file_schema = RegisterFileSchema {
            id: Uuid::new_v4().to_string(),
            run_id: run_id.map(String::from),
            sample_id: sample_id.map(String::from),
            pipeline_id: pipeline_id.map(String::from),
            description: None,
            date: now.to_string(),
            name: stored.name,
            hash: file_hash,
            fid: stored.fid,
            ftype,
            size: stored.size,
            watcher: None,
            path: stored.path,
            tier: upload_config.tier,
            retention,
            retain_until,
            legal_hold: upload_config.legal_hold,
            replicas: None,
            archived: false,
            reported_at: None,
        };

        self.api_client.register_file(file_schema)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_known_outputs() {
        let rules = CaptureRule::default_ruleset();

        assert_eq!(
            classify("sample1/sample1.cerebro.json", &rules),
            ArtefactClass::Capture { ftype: Some(FileType::CerebroModel), retention: RetentionClass::Diagnostic }
        );
        assert_eq!(
            classify("pathogen/sample1.tsv", &rules),
            ArtefactClass::Capture { ftype: Some(FileType::PathogenOutput), retention: RetentionClass::Diagnostic }
        );
        assert_eq!(
            classify("consensus/sample1.fasta", &rules),
            ArtefactClass::Capture { ftype: Some(FileType::Consensus), retention: RetentionClass::Diagnostic }
        );
        assert_eq!(
            classify("quality/fastp.json", &rules),
            ArtefactClass::Capture { ftype: Some(FileType::QualityOutput), retention: RetentionClass::Intermediate }
        );
    }

    #[test]
    fn ignores_scratch_and_logs() {
        let rules = CaptureRule::default_ruleset();
        assert_eq!(classify("work/ab/cd/scratch.tmp", &rules), ArtefactClass::Ignore);
        assert_eq!(classify("pipeline.log", &rules), ArtefactClass::Ignore);
        assert_eq!(classify(".command.sh", &rules), ArtefactClass::Ignore);
    }

    #[test]
    fn unmatched_falls_back_to_other_intermediate() {
        let rules = CaptureRule::default_ruleset();
        assert_eq!(
            classify("misc/unknown.dat", &rules),
            ArtefactClass::Capture { ftype: Some(FileType::Other), retention: RetentionClass::Intermediate }
        );
    }

    #[test]
    fn report_counts() {
        let report = CaptureReport {
            outcomes: vec![
                CaptureOutcome { relative_path: "a".into(), status: CaptureStatus::Captured { ftype: Some(FileType::Other), retention: RetentionClass::Diagnostic } },
                CaptureOutcome { relative_path: "b".into(), status: CaptureStatus::Ignored },
                CaptureOutcome { relative_path: "c".into(), status: CaptureStatus::Failed("boom".into()) },
            ],
        };
        assert_eq!(report.captured(), 1);
        assert_eq!(report.ignored(), 1);
        assert_eq!(report.failed(), 1);
        assert!(!report.ok());
    }
}
