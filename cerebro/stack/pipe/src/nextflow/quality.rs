use std::path::PathBuf;

use crate::error::WorkflowError;
use crate::parsers::nanoq::NanoqReport;
use crate::tools::scan::ScanReport;
use crate::parsers::fastp::FastpReport;

use crate::tools::umi::DeduplicationReport;
use crate::utils::{get_file_by_name, get_file_component, FileComponent};

use scrubby::report::ScrubbyReport;
use vircov::vircov::VircovSummary;



#[derive(Debug, Clone)]
pub struct QualityControlFiles {
    pub input_scan: Option<PathBuf>,
    pub reads_qc_fastp: Option<PathBuf>,
    pub reads_qc_nanoq: Option<PathBuf>,
    pub deduplication: Option<PathBuf>,
    pub internal_controls: Option<PathBuf>,
    pub synthetic_controls: Option<PathBuf>,
    pub host_depletion: Option<PathBuf>,
    pub background_alignment: Option<PathBuf>,
    pub output_scan: Option<PathBuf>,
}
impl QualityControlFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            input_scan: get_file_by_name(&path, &id, ".input.json")?,
            reads_qc_fastp: get_file_by_name(&path, &id, ".fastp.json")?,
            reads_qc_nanoq: get_file_by_name(&path, &id, ".nanoq.json")?,
            deduplication: get_file_by_name(&path, &id, ".dedup.json")?,
            internal_controls: get_file_by_name(&path, &id, ".controls.tsv")?,
            synthetic_controls: get_file_by_name(&path, &id, ".synthetic.tsv")?,
            host_depletion: get_file_by_name(&path, &id, ".host.json")?,
            background_alignment: get_file_by_name(&path, &id, ".background.tsv")?,
            output_scan: get_file_by_name(&path, &id, ".output.json")?,
        })
    }
}

pub struct QualityControlOutput {
    pub id: String,
    pub input_scan: ScanReport,
    pub reads_qc_fastp: Option<FastpReport>,
    pub reads_qc_nanoq: Option<NanoqReport>,
    pub deduplication: Option<DeduplicationReport>,
    pub internal_controls: Option<VircovSummary>,
    pub synthetic_controls: Option<VircovSummary>,
    pub host_depletion: Option<ScrubbyReport>,
    pub background_alignment: Option<VircovSummary>,
    pub output_scan: ScanReport,
}
impl QualityControlOutput {
    pub fn from(path: &PathBuf, id: Option<String>, background: bool) -> Result<Self, WorkflowError> {

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let qc_files =  QualityControlFiles::from(path, &id)?;

        if background {
            Self::get_background_qc(&id, &qc_files)
        } else {
            
            Self::get_modular_qc(&id, &qc_files)
        }
    }
    pub fn from_files(id: &str, qc_files: &QualityControlFiles, background: bool) -> Result<Self, WorkflowError> {

        if background {
            Self::get_background_qc(&id, &qc_files)
        } else {
            
            Self::get_modular_qc(&id, &qc_files)
        }
        
    }
    pub fn get_modular_qc(id: &str, qc_files: &QualityControlFiles) -> Result<Self, WorkflowError> {
        Ok(QualityControlOutput {
            id: id.to_string(), 
            input_scan: match qc_files.input_scan { 
                Some(ref path) => ScanReport::from_json(path)?, 
                None => {
                    log::error!("No required read scanning file (input) detected for {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
            reads_qc_fastp: match qc_files.reads_qc_fastp { 
                Some(ref path) => Some(FastpReport::from_json(path)?), 
                None => None
            },
            reads_qc_nanoq: match qc_files.reads_qc_nanoq { 
                Some(ref path) => Some(NanoqReport::from_json(path)?), 
                None => None
            },
            deduplication: match qc_files.deduplication { 
                Some(ref path) => Some(DeduplicationReport::from_json(path)?), 
                None => None
            },
            internal_controls: match qc_files.internal_controls { 
                Some(ref path) => Some(VircovSummary::from_tsv(path, true)?), 
                None => None
            },
            synthetic_controls: match qc_files.synthetic_controls { 
                Some(ref path) => Some(VircovSummary::from_tsv(path, true)?), 
                None => None
            },
            host_depletion: match qc_files.host_depletion { 
                Some(ref path) => Some(ScrubbyReport::from_json(path)?), 
                None => None
            },
            background_alignment:  match qc_files.background_alignment { 
                Some(ref path) => Some(VircovSummary::from_tsv(path, true)?), 
                None => None
            },
            output_scan: match qc_files.output_scan { 
                Some(ref path) => ScanReport::from_json(path)?, 
                None => {
                    log::error!("No required read scanning file (output) detected for {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
        })
    }
    pub fn get_background_qc(id: &str, qc_files: &QualityControlFiles) -> Result<Self, WorkflowError> {
        Ok(QualityControlOutput { 
            id: id.to_string(),
            input_scan: match qc_files.input_scan { 
                Some(ref path) => ScanReport::from_json(path)?, 
                None => {
                    log::error!("No required read scanning file (input) detected for {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
            reads_qc_fastp: match qc_files.reads_qc_fastp { 
                Some(ref path) => Some(FastpReport::from_json(path)?), 
                None => None
            },
            reads_qc_nanoq: match qc_files.reads_qc_nanoq { 
                Some(ref path) => Some(NanoqReport::from_json(path)?), 
                None => None
            },
            deduplication: match qc_files.deduplication { 
                Some(ref path) => Some(DeduplicationReport::from_json(path)?), 
                None => None
            },
            internal_controls: None,
            synthetic_controls: None,
            host_depletion: None,
            background_alignment: None,
            output_scan: match qc_files.output_scan { 
                Some(ref path) => ScanReport::from_json(path)?, 
                None => {
                    log::error!("No required read scanning file (output) detected for {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            }
        })
    }
}