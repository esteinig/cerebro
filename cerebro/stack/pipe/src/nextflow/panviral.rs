use std::path::PathBuf;
use crate::utils::{
    get_file_component, 
    get_file_by_name, 
    FileComponent
};

use vircov::vircov::VircovSummary;
use crate::error::WorkflowError;
use super::quality::{QualityControlFiles, QualityControlOutput};


pub struct PanviralFiles {
    pub qc: QualityControlFiles,
    pub viruses: Option<PathBuf>,
}
impl PanviralFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            qc: QualityControlFiles::from(path, id)?,
            viruses: get_file_by_name(&path, &id, ".viruses.tsv")?,
        })
    }
}

pub struct PanviralOutput {
    pub id: String,
    pub qc: QualityControlOutput,
    pub viruses: VircovSummary,
    
}
impl PanviralOutput {

    pub fn from(path: &PathBuf, id: Option<String>, background: bool) -> Result<Self, WorkflowError> {

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let files = PanviralFiles::from(&path, &id)?; 
        
        Ok(Self{
            id: id.to_string(),
            qc: QualityControlOutput::from_files(&id, &files.qc, background)? ,
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