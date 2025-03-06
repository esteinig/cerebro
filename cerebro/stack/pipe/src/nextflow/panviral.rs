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
    pub vircov: Option<PathBuf>,
}
impl PanviralFiles {
    pub fn from(path: &PathBuf, path_qc: Option<PathBuf>, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            qc: QualityControlFiles::from(match path_qc { Some(ref p) => p, None => path}, id)?,
            vircov: get_file_by_name(&path, &id, ".vircov.tsv")?,
        })
    }
}

pub struct PanviralOutput {
    pub id: String,
    pub qc: QualityControlOutput,
    pub vircov: VircovSummary,
    
}
impl PanviralOutput {

    pub fn from(path: &PathBuf, path_qc: Option<PathBuf>, id: Option<String>) -> Result<Self, WorkflowError> {

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let files = PanviralFiles::from(&path, path_qc, &id)?; 
        
        Ok(Self{
            id: id.to_string(),
            qc: QualityControlOutput::from_files(&id, &files.qc)? ,
            vircov: match files.vircov { 
                Some(ref path) => VircovSummary::from_tsv(path, true)?, 
                None => {
                    log::error!("No virus coverage file detected for: {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
        })
    }
}