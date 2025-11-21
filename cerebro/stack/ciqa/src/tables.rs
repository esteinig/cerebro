use csv::WriterBuilder;
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

fn consensus_pathogen_support_pct(
    diagnoses: &HashMap<String, Option<Diagnosis>>,
    pathogens: &HashMap<String, Option<String>>,
    consensus: &Option<String>,
) -> f64 {
    use std::collections::HashMap as Map;

    if consensus.is_none() {
        return 0.0;
    }
    let consensus_path = consensus.as_ref().unwrap();

    let mut pathogen_counts: Map<String, usize> = Map::new();
    let mut total_count = 0usize;

    for (rep_label, diag_opt) in diagnoses {
        if let Some(Diagnosis::Infectious) = diag_opt {
            if let Some(Some(path)) = pathogens.get(rep_label) {
                *pathogen_counts.entry(path.clone()).or_insert(0) += 1;
                total_count += 1;
            }
        }
    }

    if total_count == 0 {
        return 0.0;
    }

    let num = *pathogen_counts.get(consensus_path).unwrap_or(&0);
    (num as f64) * 100.0 / (total_count as f64)
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
    replicate_labels.sort_by_key(|lbl| lbl.parse::<u64>().unwrap_or(u64::MAX));

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

    // Sort sample_ids for deterministic output
    let mut sample_ids: Vec<String> = samples.keys().cloned().collect();
    sample_ids.sort();

    // diagnoses.tsv

    let mut diag_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.diagnoses)?;



    // Header: sample_id, one column per replicate, consensus, consensus_certainty
    let mut diag_header: Vec<String> = Vec::new();
    diag_header.push("sample_id".to_string());
    for rep_label in &replicate_labels {
        diag_header.push(rep_label.clone()); 
    }
    diag_header.push("consensus".to_string());
    diag_header.push("consensus_certainty".to_string()); 
    diag_writer.write_record(&diag_header)?;

    // pathogens.tsv

    let mut path_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.pathogens)?;

    let mut path_header: Vec<String> = Vec::new();
    path_header.push("sample_id".to_string());
    for rep_label in &replicate_labels {
        path_header.push(rep_label.clone()); // or format!("pathogen_{}", rep_label)
    }
    path_header.push("consensus".to_string());
    path_header.push("consensus_certainty".to_string());
    path_writer.write_record(&path_header)?;

    // rows

    for sample_id in sample_ids {
        let agg = samples.get(&sample_id).expect("sample must exist");

        // Build diagnosis row
        let mut diag_values: Vec<Option<Diagnosis>> = Vec::new();
        for rep_label in &replicate_labels {
            let value = agg.diagnoses.get(rep_label).cloned().unwrap_or(None);
            diag_values.push(value);
        }
        let diag_consensus = consensus_diagnosis(diag_values.clone());

        // Compute support percentage for consensus (0–100)
        let (mut positive_count, mut negative_count) = (0usize, 0usize);
        for d_opt in &diag_values {
            if let Some(d) = d_opt {
                match d {
                    Diagnosis::Infectious => positive_count += 1,
                    Diagnosis::NonInfectious => negative_count += 1,
                    _ => {}
                }
            }
        }
        let denom = positive_count + negative_count;
        let diag_support_pct = if let Some(ref cons) = diag_consensus {
            if denom == 0 {
                0.0
            } else {
                let num = match cons {
                    Diagnosis::Infectious => positive_count,
                    Diagnosis::NonInfectious => negative_count,
                    _ => 0,
                };
                (num as f64) * 100.0 / (denom as f64)
            }
        } else {
            0.0
        };

        let mut diag_row: Vec<String> = Vec::new();
        diag_row.push(sample_id.clone());
        for v in diag_values {
            diag_row.push(
                v.map(|d| serde_json::to_string(&d).unwrap().trim_matches('"').to_string())
                    .unwrap_or_default(),
            );
        }
        diag_row.push(
            diag_consensus
                .map(|d| serde_json::to_string(&d).unwrap().trim_matches('"').to_string())
                .unwrap_or_default(),
        );
        diag_row.push(format!("{:.2}", diag_support_pct)); 
        
        diag_writer.write_record(&diag_row)?;

        // Build pathogen row
        let mut path_row: Vec<String> = Vec::new();
        path_row.push(sample_id.clone());
        for rep_label in &replicate_labels {
            let value = agg.pathogens.get(rep_label).cloned().unwrap_or(None);
            path_row.push(value.unwrap_or_default());
        }

        let pathogen_consensus = consensus_pathogen(&agg.diagnoses, &agg.pathogens);
        let pathogen_support_pct =  consensus_pathogen_support_pct(&agg.diagnoses, &agg.pathogens, &pathogen_consensus);

        path_row.push(pathogen_consensus.unwrap_or_default());
        path_row.push(format!("{:.2}", pathogen_support_pct));

        path_writer.write_record(&path_row)?;
    }

    diag_writer.flush()?;
    path_writer.flush()?;

    Ok(())
}