//! Run finalisation orchestration.
//!
//! When a pipeline run completes, its outputs must be captured into cerebro-fs —
//! and a provenance manifest sealed and stored — **before** the Nextflow
//! execution directory is cleaned up, so nothing is lost to scratch reclamation.
//! This module performs that capture using the tower's [`FileSystemClient`] hook
//! (FS-5), reusing the capture and manifest building.
//!
//! Capture targets are chosen per layout: if the execution directory contains a
//! per-sample subdirectory (named exactly the sample id), each sample's outputs
//! are captured and linked to that sample; otherwise the whole execution
//! directory is captured once under the run. The classifier already ignores the
//! Nextflow `work/` scratch, so pointing it at the execution directory is safe.
//!
//! The blocking cerebro-fs work runs on a dedicated blocking thread; the caller
//! gates cleanup on [`FinalizeReport::ok`].

use std::path::{Path, PathBuf};

use cerebro_fs::capture::CaptureRule;
use cerebro_fs::client::{FileSystemClient, UploadConfig};
use cerebro_model::api::files::manifest::{ManifestProvenance, ReferenceDb, ToolVersion};
use cerebro_model::api::files::retention::RetentionPolicy;
use cerebro_model::api::stage::model::StagedSample;

use crate::error::TowerError;

/// Summary of a run finalisation.
#[derive(Debug, Clone)]
pub struct FinalizeReport {
    /// Number of distinct runs finalised.
    pub runs: usize,
    /// Number of capture targets (per-sample or per-run) successfully captured.
    pub captured: usize,
    /// Human-readable failures, one per target that did not fully capture.
    pub failures: Vec<String>,
}
impl FinalizeReport {
    /// True when every target captured cleanly — the precondition for cleanup.
    pub fn ok(&self) -> bool {
        self.failures.is_empty()
    }
}

/// Distinct values in first-seen order.
fn distinct_in_order(values: impl Iterator<Item = String>) -> Vec<String> {
    let mut seen: Vec<String> = Vec::new();
    for value in values {
        if !seen.iter().any(|v| v == &value) {
            seen.push(value);
        }
    }
    seen
}

/// Distinct run identifiers in first-seen order.
fn distinct_run_ids(samples: &[StagedSample]) -> Vec<String> {
    distinct_in_order(samples.iter().map(|s| s.run_id.clone()))
}

/// Candidate directories that may hold pipeline-emitted metadata for `dir`: the
/// dir itself, its `output/` child, and the same two for the parent (so per-sample
/// capture still finds run-level sidecars).
fn provenance_dirs(dir: &Path) -> Vec<PathBuf> {
    let mut dirs = vec![dir.to_path_buf(), dir.join("output")];
    if let Some(parent) = dir.parent() {
        dirs.push(parent.to_path_buf());
        dirs.push(parent.join("output"));
    }
    dirs
}

/// Load and assemble pipeline-emitted run provenance.
///
/// Base metadata comes from `provenance.json` (pipeline name/version/params).
/// `tool_versions` are assembled from `tool_versions.tsv` (`name<TAB>version`) and
/// `reference_dbs` from `reference_dbs.tsv` (`name<TAB>version[<TAB>hash]`) when the
/// JSON didn't already carry them, so the sealed manifest records the real tools
/// and database builds used. Missing/invalid files default to empty.
fn load_provenance(dir: &Path) -> ManifestProvenance {
    let dirs = provenance_dirs(dir);

    let mut provenance = dirs
        .iter()
        .find_map(|d| {
            std::fs::read(d.join("provenance.json"))
                .ok()
                .and_then(|bytes| serde_json::from_slice::<ManifestProvenance>(&bytes).ok())
        })
        .unwrap_or_default();

    if provenance.tool_versions.is_empty() {
        if let Some(versions) = dirs
            .iter()
            .find_map(|d| read_tool_versions(&d.join("tool_versions.tsv")))
        {
            provenance.tool_versions = versions;
        }
    }
    if provenance.reference_dbs.is_empty() {
        if let Some(dbs) = dirs
            .iter()
            .find_map(|d| read_reference_dbs(&d.join("reference_dbs.tsv")))
        {
            provenance.reference_dbs = dbs;
        }
    }

    provenance
}

/// Parse a `name<TAB>version` TSV into [`ToolVersion`] entries (`None` if absent
/// or empty).
fn read_tool_versions(path: &Path) -> Option<Vec<ToolVersion>> {
    let text = std::fs::read_to_string(path).ok()?;
    let versions: Vec<ToolVersion> = text
        .lines()
        .filter_map(|line| {
            let mut it = line.splitn(2, '\t');
            let name = it.next()?.trim();
            let version = it.next().unwrap_or("").trim();
            if name.is_empty() {
                None
            } else {
                Some(ToolVersion {
                    name: name.to_string(),
                    version: version.to_string(),
                })
            }
        })
        .collect();
    if versions.is_empty() {
        None
    } else {
        Some(versions)
    }
}

/// Parse a `name<TAB>version[<TAB>hash]` TSV into [`ReferenceDb`] entries (`None`
/// if absent or empty).
fn read_reference_dbs(path: &Path) -> Option<Vec<ReferenceDb>> {
    let text = std::fs::read_to_string(path).ok()?;
    let dbs: Vec<ReferenceDb> = text
        .lines()
        .filter_map(|line| {
            let mut it = line.split('\t');
            let name = it.next()?.trim();
            let version = it.next().unwrap_or("").trim();
            let hash = it
                .next()
                .map(|h| h.trim().to_string())
                .filter(|h| !h.is_empty());
            if name.is_empty() {
                None
            } else {
                Some(ReferenceDb {
                    name: name.to_string(),
                    version: version.to_string(),
                    hash,
                })
            }
        })
        .collect();
    if dbs.is_empty() {
        None
    } else {
        Some(dbs)
    }
}

/// Capture a completed run's outputs and provenance into cerebro-fs, then return
/// a report the caller uses to gate execution-directory cleanup.
///
/// Blocking cerebro-fs operations run on a blocking thread so the async runtime
/// is not stalled.
pub async fn finalize_run(
    fs_client: FileSystemClient,
    execution_dir: PathBuf,
    staged_samples: Vec<StagedSample>,
) -> Result<FinalizeReport, TowerError> {
    tokio::task::spawn_blocking(move || {
        finalize_run_blocking(&fs_client, &execution_dir, &staged_samples)
    })
    .await
    .map_err(|e| TowerError::Finalize(format!("capture task panicked: {e}")))?
}

fn finalize_run_blocking(
    fs_client: &FileSystemClient,
    execution_dir: &Path,
    staged_samples: &[StagedSample],
) -> Result<FinalizeReport, TowerError> {
    let upload_config = UploadConfig {
        retention_policy: RetentionPolicy::from_env(),
        ..UploadConfig::default()
    };
    let rules = CaptureRule::default_ruleset();

    let run_ids = distinct_run_ids(staged_samples);
    let mut captured = 0usize;
    let mut failures = Vec::new();

    for run_id in &run_ids {
        let run_samples: Vec<&StagedSample> = staged_samples
            .iter()
            .filter(|s| &s.run_id == run_id)
            .collect();

        // Per-sample capture where a per-sample output subdirectory exists.
        let mut captured_per_sample = false;
        for sample in &run_samples {
            let sample_dir = execution_dir.join(&sample.sample_id);
            if sample_dir.is_dir() {
                match capture_target(
                    fs_client,
                    &sample_dir,
                    run_id,
                    Some(&sample.sample_id),
                    &rules,
                    &upload_config,
                ) {
                    Ok(()) => {
                        captured += 1;
                        captured_per_sample = true;
                    }
                    Err(e) => failures.push(format!("{}/{}: {}", run_id, sample.sample_id, e)),
                }
            }
        }

        // Otherwise capture the whole execution directory once under the run.
        if !captured_per_sample {
            match capture_target(
                fs_client,
                execution_dir,
                run_id,
                None,
                &rules,
                &upload_config,
            ) {
                Ok(()) => captured += 1,
                Err(e) => failures.push(format!("{}: {}", run_id, e)),
            }
        }
    }

    Ok(FinalizeReport {
        runs: run_ids.len(),
        captured,
        failures,
    })
}

/// Capture one target directory and seal+store its provenance manifest.
fn capture_target(
    fs_client: &FileSystemClient,
    dir: &Path,
    run_id: &str,
    sample_id: Option<&str>,
    rules: &[CaptureRule],
    upload_config: &UploadConfig,
) -> Result<(), TowerError> {
    let report = fs_client.capture_outputs(
        Some(run_id.to_string()),
        sample_id.map(String::from),
        None,
        dir,
        rules,
        upload_config,
    )?;
    if !report.ok() {
        return Err(TowerError::Finalize(format!(
            "{} of {} file(s) failed to capture",
            report.failed(),
            report.outcomes.len()
        )));
    }

    let provenance = load_provenance(dir);
    let manifest = fs_client.build_run_manifest(
        Some(run_id.to_string()),
        sample_id.map(String::from),
        None,
        provenance,
    )?;
    fs_client.capture_manifest(&manifest, upload_config)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distinct_in_order_dedupes_first_seen() {
        let ids = vec!["RUN2".to_string(), "RUN1".to_string(), "RUN2".to_string()];
        assert_eq!(
            distinct_in_order(ids.into_iter()),
            vec!["RUN2".to_string(), "RUN1".to_string()]
        );
    }

    #[test]
    fn report_ok_only_without_failures() {
        let ok = FinalizeReport {
            runs: 1,
            captured: 1,
            failures: vec![],
        };
        assert!(ok.ok());
        let bad = FinalizeReport {
            runs: 1,
            captured: 0,
            failures: vec!["x".into()],
        };
        assert!(!bad.ok());
    }
}
