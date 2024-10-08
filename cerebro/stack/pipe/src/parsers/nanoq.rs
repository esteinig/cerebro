
use std::{fs::File, io::BufReader, path::PathBuf};

use serde::{Deserialize, Serialize};

use std::io::Write;
use crate::error::WorkflowError;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NanoqReport {
    pub reads: u64,
    pub bases: u64,
    pub n50: u64,
    pub longest: u64,
    pub shortest: u64,
    pub mean_length: u64,
    pub median_length: u64,
    pub mean_quality: Option<f64>,
    pub median_quality: Option<f64>,
    pub filtered: u64
}

impl NanoqReport {
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