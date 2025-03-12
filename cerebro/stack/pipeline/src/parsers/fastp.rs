use std::{fs::File, io::BufReader, path::PathBuf};

use serde::{Deserialize, Serialize};

use std::io::Write;
use crate::error::WorkflowError;


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastpReport {
    pub summary: Summary,
    #[serde(rename = "filtering_result")]
    pub filter: Filter,
    // May not be run with duplication option
    pub duplication: Option<Duplication>,
    #[serde(rename = "adapter_cutting")]
    pub adapter: Option<Adapters>
}

impl FastpReport {
    pub fn to_json(report: &Self, output: &PathBuf) -> Result<(), WorkflowError> {
        let mut file = std::fs::File::create(output)?;
        let json_string = serde_json::to_string_pretty(report)?;
        file.write_all(json_string.as_bytes())?;
        Ok(())
    }
    pub fn from_json(report: &PathBuf) -> Result<Self, WorkflowError> {
        let mut reader = BufReader::new(File::open(&report)?);
        let report: Self = serde_json::from_reader(&mut reader)?;
        Ok(report)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Summary {
    #[serde(rename = "fastp_version")] 
    pub version: String,
    #[serde(rename = "before_filtering")] 
    pub before: ReadSummary,
    #[serde(rename = "after_filtering")] 
    pub after: ReadSummary
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadSummary {
    #[serde(rename = "total_reads")] 
    pub reads: u64,
    #[serde(rename = "total_bases")] 
    pub bases: u64,
    #[serde(rename = "q20_rate")] 
    pub q20: f64,
    #[serde(rename = "q30_rate")] 
    pub q30: f64,
    #[serde(rename = "read1_mean_length")] 
    pub mean_length_r1: u64,
    #[serde(rename = "read2_mean_length")] 
    pub mean_length_r2: u64,
    #[serde(rename = "gc_content")] 
    pub gc: f64
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Filter {
    #[serde(rename = "passed_filter_reads")] 
    pub pass_filter: u64,
    #[serde(rename = "low_quality_reads")] 
    pub low_quality: u64,
    #[serde(rename = "low_complexity_reads")] 
    pub low_complexity:  Option<u64>,
    #[serde(rename = "too_many_N_reads")] 
    pub min_missing: u64,
    #[serde(rename = "too_short_reads")] 
    pub min_length: u64,
    #[serde(rename = "too_long_reads")] 
    pub max_length: u64
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Duplication {
    pub rate: f64
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Adapters {
    #[serde(rename = "adapter_trimmed_reads")] 
    pub total_reads: u64,
    #[serde(rename = "adapter_trimmed_bases")] 
    pub total_bases: u64,
    #[serde(rename = "read1_adapter_sequence")] 
    pub sequence_r1: String,
    #[serde(rename = "read2_adapter_sequence")] 
    pub sequence_r2: String
}

