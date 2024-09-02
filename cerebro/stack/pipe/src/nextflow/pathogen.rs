use std::path::PathBuf;

use crate::error::WorkflowError;
use crate::tools::scan::ScanReport;
use crate::parsers::fastp::FastpReport;

use crate::tools::umi::DeduplicationReport;
use crate::utils::{
    get_file_component, 
    get_file_by_name, 
    FileComponent
};

use scrubby::report::ScrubbyReport;
use vircov::vircov::VircovSummary;

#[derive(Debug, Clone)]
pub struct PathogenFiles {
    pub qc: QualityControlFiles,

}
impl PathogenFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            qc: QualityControlFiles::from(path, id)?
        })
    }
}

#[derive(Debug, Clone)]
pub struct QualityControlFiles {
    pub input_scan: Option<PathBuf>,
    pub reads_qc: Option<PathBuf>,
    pub deduplication: Option<PathBuf>,
    pub internal_controls: Option<PathBuf>,
    pub synthetic_controls: Option<PathBuf>,
    pub host_depletion: Option<PathBuf>,
    pub background_depletion: Option<PathBuf>,
    pub output_scan: Option<PathBuf>,
}
impl QualityControlFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            input_scan: get_file_by_name(&path, &id, ".input.json")?,
            reads_qc: get_file_by_name(&path, &id, ".reads.json")?,
            deduplication: get_file_by_name(&path, &id, ".dedup.json")?,
            internal_controls: get_file_by_name(&path, &id, ".controls.tsv")?,
            synthetic_controls: get_file_by_name(&path, &id, ".synthetic.tsv")?,
            host_depletion: get_file_by_name(&path, &id, ".host.json")?,
            background_depletion: get_file_by_name(&path, &id, ".background.tsv")?,
            output_scan: get_file_by_name(&path, &id, ".output.json")?,
        })
    }
}

pub struct PathogenOutput {
    pub id: String,
    pub qc: QualityControlOutput,
    
}
impl PathogenOutput {

    pub fn from(path: &PathBuf, id: Option<String>) -> Result<Self, WorkflowError> {

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let files = QualityControlFiles::from(&path, &id)?;

        log::info!("Quality control files: {:#?}", files);
        
        Ok(Self{
            id: id.to_string(),
            qc: QualityControlOutput { 
                input_scan: match files.input_scan { 
                    Some(ref path) => ScanReport::from_json(path)?, 
                    None => {
                        log::error!("No required read scanning file (input) detected for {id}");
                        return Err(WorkflowError::PipelineOutputNotFound)
                    }
                },
                reads_qc: match files.reads_qc { 
                    Some(ref path) => Some(FastpReport::from_json(path)?), 
                    None => None
                },
                deduplication: match files.deduplication { 
                    Some(ref path) => Some(DeduplicationReport::from_json(path)?), 
                    None => None
                },
                internal_controls: match files.internal_controls { 
                    Some(ref path) => Some(VircovSummary::from_tsv(path, true)?), 
                    None => None
                },
                synthetic_controls: match files.synthetic_controls { 
                    Some(ref path) => Some(VircovSummary::from_tsv(path, true)?), 
                    None => None
                },
                host_depletion: match files.host_depletion { 
                    Some(ref path) => Some(ScrubbyReport::from_json(path)?), 
                    None => None
                },
                background_depletion:  match files.background_depletion { 
                    Some(ref path) => Some(VircovSummary::from_tsv(path, true)?), 
                    None => None
                },
                output_scan: match files.output_scan { 
                    Some(ref path) => ScanReport::from_json(path)?, 
                    None => {
                        log::error!("No required read scanning file (output) detected for {id}");
                        return Err(WorkflowError::PipelineOutputNotFound)
                    }
                },
            }
        })
    }
}

pub struct QualityControlOutput {
    pub input_scan: ScanReport,
    pub reads_qc: Option<FastpReport>,
    pub deduplication: Option<DeduplicationReport>,
    pub internal_controls: Option<VircovSummary>,
    pub synthetic_controls: Option<VircovSummary>,
    pub host_depletion: Option<ScrubbyReport>,
    pub background_depletion: Option<VircovSummary>,
    pub output_scan: ScanReport,
}
