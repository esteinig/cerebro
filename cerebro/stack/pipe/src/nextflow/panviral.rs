use std::path::PathBuf;

use crate::{error::WorkflowError, parsers::fastp::FastpReport, utils::{get_file_component, FileComponent}};

use scrubby::report::ScrubbyReport;
use vircov::vircov::VircovSummary;

use super::utils::get_file_by_name;

pub struct PanviralFiles {
    pub reads: Option<PathBuf>,
    pub controls: Option<PathBuf>,
    pub host: Option<PathBuf>,
    pub viruses: Option<PathBuf>,
}
impl PanviralFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            reads: get_file_by_name(&path, &id, ".reads.json")?,
            host: get_file_by_name(&path, &id, ".host.json")?,
            viruses: get_file_by_name(&path, &id, ".viruses.tsv")?,
            controls: get_file_by_name(&path, &id, ".controls.tsv")?,
        })
    }
}

pub struct PanviralOutput {
    pub id: String,
    pub reads: FastpReport,
    pub controls: Option<VircovSummary>,
    pub host: Option<ScrubbyReport>,
    pub viruses: VircovSummary,
    
}
impl PanviralOutput {

    pub fn from(path: &PathBuf, id: Option<String>) -> Result<Self, WorkflowError> {

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let files = PanviralFiles::from(&path, &id)?;
        
        Ok(Self{
            id: id.to_string(),
            reads: match files.reads { 
                Some(ref path) => FastpReport::from_json(path)?, 
                None => {
                    log::error!("No read quality control file detected for {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
            host: match files.host { 
                Some(ref path) => Some(ScrubbyReport::from_json(path)?), 
                None => None
            },
            controls: match files.controls { 
                Some(ref path) => Some(VircovSummary::from_tsv(path, true)?), 
                None => None
            },
            viruses: match files.viruses { 
                Some(ref path) => VircovSummary::from_tsv(path, true)?, 
                None => {
                    log::error!("No virus coverage file detected for {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
        })
    }
}