use meta_gpt::gpt::{Diagnosis, DiagnosticResult};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, BTreeMap};
use std::fs;
use std::path::{Path, PathBuf};

use crate::terminal::PredictionSummaryArgs;
use crate::utils::write_tsv;

// One row per sample for diagnoses.tsv
#[derive(Debug, Serialize)]
struct DiagnosisRow {
    pub sample_id: String,
    #[serde(flatten)]
    pub replicates: BTreeMap<String, Option<Diagnosis>>,
    pub consensus: Option<Diagnosis>,
}

// One row per sample for pathogens.tsv
#[derive(Debug, Serialize)]
struct PathogenRow {
    pub sample_id: String,
    #[serde(flatten)]
    pub replicates: BTreeMap<String, Option<String>>,
    pub consensus: Option<String>,
}

// Internal aggregation for each sample
#[derive(Debug, Default)]
struct SampleAggregate {
    diagnoses: HashMap<String, Option<Diagnosis>>, // replicate label -> diagnosis
    pathogens: HashMap<String, Option<String>>,    // replicate label -> pathogen
}

fn replicate_label_for_dir(dir: &Path, fallback_index: usize) -> String {
    let name = dir
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or_else(|| panic!("Invalid UTF-8 in directory name: {:?}", dir));

    if let Some((_, suffix)) = name.rsplit_once('_') {
        if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
            return suffix.to_string();
        }
    }

    // Fallback: order of the folders, 1-based
    fallback_index.to_string()
}

fn sample_id_from_path(path: &Path) -> Option<String> {
    let stem = path.file_stem()?.to_str()?;
    // Pattern {sample_id}.*.json - take part before first '.'
    let sample_id = stem.split_once('.').map(|(id, _)| id).unwrap_or(stem);
    Some(sample_id.to_string())
}

fn consensus_diagnosis(diags: Vec<Option<Diagnosis>>) -> Option<Diagnosis> {
    let mut positive_count = 0;
    let mut negative_count = 0;

    for d_opt in diags {
        if let Some(d) = d_opt {
            match d {
                Diagnosis::Infectious => positive_count += 1,
                Diagnosis::NonInfectious => negative_count += 1,
                // Everything else ignored
                Diagnosis::InfectiousReview
                | Diagnosis::NonInfectiousReview
                | Diagnosis::Tumor
                | Diagnosis::Unknown => {}
            }
        }
    }

    // Mimic SampleReview logic: tie or positives > negatives → Positive
    if positive_count >= negative_count {
        Some(Diagnosis::Infectious)
    } else {
        Some(Diagnosis::NonInfectious)
    }
}

fn consensus_pathogen(
    diagnoses: &HashMap<String, Option<Diagnosis>>,
    pathogens: &HashMap<String, Option<String>>,
) -> Option<String> {
    use std::collections::HashMap as Map;

    let mut positive_count = 0;
    let mut negative_count = 0;
    let mut pathogen_counts: Map<String, usize> = Map::new();

    for (rep_label, diag_opt) in diagnoses {
        if let Some(d) = diag_opt {
            match d {
                Diagnosis::Infectious => {
                    positive_count += 1;
                    if let Some(Some(path)) = pathogens.get(rep_label) {
                        *pathogen_counts.entry(path.clone()).or_insert(0) += 1;
                    }
                }
                Diagnosis::NonInfectious => {
                    negative_count += 1;
                }
                Diagnosis::InfectiousReview
                | Diagnosis::NonInfectiousReview
                | Diagnosis::Tumor
                | Diagnosis::Unknown => {}
            }
        }
    }

    let consensus_diag = if positive_count >= negative_count {
        Some(Diagnosis::Infectious)
    } else {
        Some(Diagnosis::NonInfectious)
    };

    if consensus_diag == Some(Diagnosis::Infectious) {
        pathogen_counts
            .into_iter()
            .max_by_key(|(_, count)| *count)
            .map(|(pathogen, _)| pathogen)
    } else {
        None
    }
}



pub fn summarize_predictions(args: &PredictionSummaryArgs) -> anyhow::Result<()> {

    // Determine replicate labels, preserving order of input_dirs
    let mut replicate_labels: Vec<String> = Vec::new();
    let mut replicate_dirs: Vec<PathBuf> = Vec::new();

    for (i, dir) in args.input_dirs.iter().enumerate() {
        let label = replicate_label_for_dir(dir, i + 1);
        replicate_labels.push(label);
        replicate_dirs.push(dir.to_path_buf());
    }

    // Aggregate all DiagnosticResult per {sample_id} across replicates
    let mut samples: HashMap<String, SampleAggregate> = HashMap::new();

    for (rep_idx, dir) in replicate_dirs.iter().enumerate() {
        let rep_label = &replicate_labels[rep_idx];

        for entry in fs::read_dir(dir)? {
            
            let entry = entry?;
            let path = entry.path();

            if path.extension().and_then(|s| s.to_str()) != Some("json") {
                continue;
            }

            let Some(sample_id) = sample_id_from_path(&path) else {
                continue;
            };

            let file = fs::File::open(&path)?;
            let result: DiagnosticResult = serde_json::from_reader(file)?;

            let agg = samples
                .entry(sample_id)
                .or_insert_with(SampleAggregate::default);

            // Diagnosis is always present in DiagnosticResult.
            agg.diagnoses.insert(rep_label.clone(), Some(result.diagnosis));
            // Pathogen is Option<String>.
            agg.pathogens.insert(rep_label.clone(), result.pathogen);
        }
    }

    // Build diagnosis rows
    let mut diagnosis_rows: Vec<DiagnosisRow> = Vec::new();
    let mut pathogen_rows: Vec<PathogenRow> = Vec::new();

    // Sort sample_ids for deterministic output
    let mut sample_ids: Vec<String> = samples.keys().cloned().collect();
    sample_ids.sort();

    for sample_id in sample_ids {
        let agg = samples.get(&sample_id).expect("sample must exist");

        // Build per-replicate columns for diagnoses
        let mut diag_map: BTreeMap<String, Option<Diagnosis>> = BTreeMap::new();
        for rep_label in &replicate_labels {
            let value = agg.diagnoses.get(rep_label).cloned().unwrap_or(None);
            diag_map.insert(rep_label.clone(), value);
        }

        let mut diag = Vec::new();
        for i in diag_map.iter() {
            diag.push(i.1.clone())
        }

        let diag_consensus = consensus_diagnosis(diag);

        diagnosis_rows.push(DiagnosisRow {
            sample_id: sample_id.clone(),
            replicates: diag_map,
            consensus: diag_consensus,
        });

        // Build per-replicate columns for pathogens
        let mut path_map: BTreeMap<String, Option<String>> = BTreeMap::new();
        for rep_label in &replicate_labels {
            let value = agg.pathogens.get(rep_label).cloned().unwrap_or(None);
            path_map.insert(rep_label.clone(), value);
        }

        let pathogen_consensus = consensus_pathogen(&agg.diagnoses, &agg.pathogens);

        pathogen_rows.push(PathogenRow {
            sample_id,
            replicates: path_map,
            consensus: pathogen_consensus,
        });
    }

    // Write TSVs using your existing helper
    write_tsv(&diagnosis_rows, &args.diagnoses, true)?;
    write_tsv(&pathogen_rows, &args.pathogens, true)?;

    Ok(())
}