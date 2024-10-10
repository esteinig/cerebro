use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}, path::Path};

use serde::{Deserialize, Serialize};

use crate::{error::WorkflowError, nextflow::pathogen::PathogenOutput, utils::{read_tsv, write_tsv}};

use super::quality::QualityControl;

fn compute_rpm(reads: u64, input_reads: u64) -> Option<f64> {
    if input_reads == 0 {
        None  // Cannot compute RPM when input_reads is 0
    } else {
        Some((reads as f64 / input_reads as f64) * 1_000_000.0)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetectionFilter {
    pub taxids: Option<Vec<String>>,
    pub names: Option<Vec<String>>,
    pub ranks: Option<Vec<String>>
}
impl PathogenDetectionFilter {
    pub fn new(
        taxids: Option<Vec<String>>,
        names: Option<Vec<String>>,
        ranks: Option<Vec<String>>
    ) -> Result<Self, WorkflowError> {
        Ok(Self {
            taxids, names, ranks
        })
    }
    pub fn from_json(path: &Path) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(path)?);
        let pathogen = serde_json::from_reader(reader)?;
        Ok(pathogen)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetection {
    pub id: String,
    pub records: Vec<PathogenDetectionRecord>
}
impl PathogenDetection {
    pub fn from_pathogen(output: &PathogenOutput, quality: &QualityControl) -> Self {

        let input_reads = quality.reads.input_reads;
        let mut detection_map: HashMap<String, PathogenDetectionRecord> = HashMap::new();

        // Process KrakenReport
        if let Some(kraken_report) = &output.profile.kraken2 {
            for record in &kraken_report.records {

                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                    .or_insert(PathogenDetectionRecord {
                        taxid,
                        kraken_reads: Some(record.reads),
                        kraken_rpm: compute_rpm(record.reads, input_reads),
                        kraken_name: Some(record.taxname.trim().to_string()),
                        kraken_rank: Some(record.tax_level.clone()),
                        metabuli_reads: None,
                        metabuli_rpm: None,
                        metabuli_name: None,
                        metabuli_rank: None,
                        bracken_reads: None,
                        bracken_rpm: None,
                        bracken_name: None,
                        bracken_rank: None
                    });
                entry.kraken_reads = Some(record.reads);
                entry.kraken_name = Some(record.taxname.trim().to_string());
                entry.kraken_rpm = compute_rpm(record.reads_direct, input_reads);
                entry.kraken_rank = Some(record.tax_level.clone())
            }
        }

        // Process MetabuliReport
        if let Some(metabuli_report) = &output.profile.metabuli {
            for record in &metabuli_report.records {

                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                    .or_insert(PathogenDetectionRecord {
                        taxid,
                        kraken_reads: None,
                        kraken_rank: None,
                        kraken_name: None,
                        kraken_rpm: None,
                        metabuli_reads: Some(record.reads),
                        metabuli_rpm: compute_rpm(record.reads, input_reads),
                        metabuli_name: Some(record.taxname.trim().to_string()),
                        metabuli_rank: Some(record.tax_level.clone()),
                        bracken_reads: None,
                        bracken_rpm: None,
                        bracken_name: None,
                        bracken_rank: None
                        
                    });
                entry.metabuli_reads = Some(record.reads);
                entry.metabuli_name = Some(record.taxname.trim().to_string());
                entry.metabuli_rpm = compute_rpm(record.reads_direct, input_reads);
                entry.metabuli_rank = Some(record.tax_level.clone());
            }
        }

        // Process BrackenReport
        if let Some(bracken_report) = &output.profile.bracken {
            for record in &bracken_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                    .or_insert(PathogenDetectionRecord {
                        taxid,
                        kraken_reads: None,
                        kraken_rank: None,
                        kraken_name: None,
                        kraken_rpm: None,
                        metabuli_reads: None,
                        metabuli_rpm: None,
                        metabuli_name: None,
                        metabuli_rank: None,
                        bracken_reads: Some(record.reads),
                        bracken_rpm: compute_rpm(record.reads, input_reads),
                        bracken_name: Some(record.taxname.trim().to_string()),
                        bracken_rank: Some(record.tax_level.clone())
                    });
                entry.bracken_reads = Some(record.reads);
                entry.bracken_name = Some(record.taxname.trim().to_string());
                entry.bracken_rpm = compute_rpm(record.reads, input_reads);
                entry.bracken_rank = Some(record.tax_level.clone());
            }
        }

        let records = detection_map
            .into_iter()
            .map(|(_, record)| record)
            .collect();

        Self { id: output.id.to_string(), records }
    }
    pub fn filter_by_taxonomy(
        &self,
        taxids: Option<Vec<String>>, 
        names: Option<Vec<String>>, 
        ranks: Option<Vec<String>>
    ) -> Vec<PathogenDetectionRecord> {
        // Filter records by taxid, name, or rank
        self.records.iter().filter(|record| {
            // Check if the taxid matches
            let taxid_match = if let Some(taxids) = &taxids {
                taxids.contains(&record.taxid)
            } else {
                true
            };

            // Check if any name matches (from kraken, bracken, or metabuli)
            let name_match = if let Some(names) = &names {
                let kraken_name_match = record.kraken_name.as_ref().map_or(false, |name| names.contains(name));
                let bracken_name_match = record.bracken_name.as_ref().map_or(false, |name| names.contains(name));
                let metabuli_name_match = record.metabuli_name.as_ref().map_or(false, |name| names.contains(name));
                kraken_name_match || bracken_name_match || metabuli_name_match
            } else {
                true
            };

            // Check if any rank matches (from kraken, bracken, or metabuli)
            let rank_match = if let Some(ranks) = &ranks {
                let kraken_rank_match = record.kraken_rank.as_ref().map_or(false, |rank| ranks.contains(rank));
                let bracken_rank_match = record.bracken_rank.as_ref().map_or(false, |rank| ranks.contains(rank));
                let metabuli_rank_match = record.metabuli_rank.as_ref().map_or(false, |rank| ranks.contains(rank));
                kraken_rank_match || bracken_rank_match || metabuli_rank_match
            } else {
                true
            };

            // Record should be included if it matches taxid, name, or rank criteria
            taxid_match && name_match && rank_match
        })
        .cloned()
        .collect()
    }

    pub fn write_records(records: &Vec<PathogenDetectionRecord>, path: &Path) -> Result<(), WorkflowError> {
        write_tsv(&records, path, true)
    }
    pub fn read_records(path: &Path) -> Result<Vec<PathogenDetectionRecord>, WorkflowError> {
        read_tsv(path, false, true)
    }
    pub fn to_tsv(&self, path: &Path) -> Result<(), WorkflowError> {
        write_tsv(&self.records, path, true)
    }
    pub fn from_tsv(&self, path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self {
            id: id.to_string(),
            records: read_tsv(path, false, true)?
        })
    }
    pub fn to_json(&self, path: &Path) -> Result<(), WorkflowError> {
        let writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }
    pub fn from_json(path: &Path) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(path)?);
        let pathogen = serde_json::from_reader(reader)?;
        Ok(pathogen)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetectionRecord {
    // pub id: String,
    // pub rank: String,
    pub taxid: String,
    // pub input_reads: u64,
    // pub qc_reads: u64,
    pub kraken_name: Option<String>,
    pub kraken_rank: Option<String>,
    pub kraken_reads: Option<u64>,
    pub kraken_rpm: Option<f64>,
    pub bracken_name: Option<String>,
    pub bracken_rank: Option<String>,
    pub bracken_reads: Option<u64>,
    pub bracken_rpm: Option<f64>,
    pub metabuli_name: Option<String>,
    pub metabuli_rank: Option<String>,
    pub metabuli_reads: Option<u64>,
    pub metabuli_rpm: Option<f64>,
    // pub kmcp_reads: Option<u64>,
    // pub scan_reads: Option<u64>,
    // pub remap_reads: Option<u64>
}