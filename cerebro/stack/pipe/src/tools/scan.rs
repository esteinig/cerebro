
use std::io::Write;
use rayon::iter::{ParallelIterator, IntoParallelRefIterator};
use serde::{Deserialize, Serialize};
use std::{fs::File, io::BufReader, path::PathBuf};

use crate::{error::WorkflowError, utils::parse_fastx_file_with_check};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScanReport {
    pub reads: usize,
    pub bases: usize
}
impl ScanReport {
    pub fn new(reads: usize, bases: usize) -> Self {
        Self { reads, bases }
    }
    pub fn from_vec(reports: &Vec<Self>) -> Self {
        let mut combined = Self::new(0, 0);
        for report in reports {
            combined.reads += report.reads;
            combined.bases += report.bases;
        }
        combined
    }
    pub fn from_json(path: &PathBuf) -> Result<Self, WorkflowError> {
        let mut reader = BufReader::new(File::open(&path)?);
        let report: Self = serde_json::from_reader(&mut reader)?;
        Ok(report)
    }
    pub fn to_json(&self, path: &PathBuf) -> Result<(), WorkflowError> {
        let mut file = std::fs::File::create(path)?;
        let json_string = serde_json::to_string_pretty(self)?;
        file.write_all(json_string.as_bytes())?;
        Ok(())
    }
}

pub struct ScanReads {
    pub input: Vec<PathBuf>,
    pub report: ScanReport,
}
impl ScanReads {
    pub fn new(input: &Vec<PathBuf>) -> Self {
        Self {
            input: input.to_owned(),
            report: ScanReport::new(0, 0)
        }
    }
    pub fn report(&self) -> Result<ScanReport, WorkflowError> {

        let reports = self.input.par_iter().map(|path| {

            let mut reads = 0;
            let mut bases = 0;

           match parse_fastx_file_with_check(&path) {
                Ok(reader) => {
                    if let Some(mut reader) = reader {
                        while let Some(record) = reader.next() {
                            match record {
                                Ok(rec) => {
                                    reads += 1;
                                    bases += rec.num_bases();
                                },
                                Err(err) => {
                                    log::error!("{}", err.to_string())
                                }
                            };
                        }
                    } else {
                        log::warn!(
                            "Input file is empty: {}", 
                            path.display()
                        )
                    }
                },
                Err(err) => {
                    log::error!("{}", err.to_string())
                }
            };
            
            ScanReport::new(reads, bases)
        }).collect::<Vec<ScanReport>>();

        Ok(ScanReport::from_vec(&reports))
    }
}

