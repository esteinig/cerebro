use std::path::PathBuf;
use crate::utils::{
    get_file_component, 
    get_file_by_name, 
    FileComponent
};

use vircov::vircov::VircovSummary;
use crate::error::WorkflowError;
use super::quality::{QualityControlFiles, QualityControlOutput};


pub struct MagFiles {
    pub qc: QualityControlFiles,
    pub ncbi_blast: Option<PathBuf>,
}
impl MagFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            qc: QualityControlFiles::from(path, id)?,
            ncbi_blast: get_file_by_name(&path, &id, ".vircov.tsv")?,
        })
    }
}

pub struct MagOutput {
    pub id: String,
    pub qc: QualityControlOutput,
    pub ncbi_blast: String,
    
}
impl MagOutput {

    pub fn from(path: &PathBuf, id: Option<String>) -> Result<Self, WorkflowError> {

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let files = MagFiles::from(&path, &id)?; 
        
        Ok(Self{
            id: id.to_string(),
            qc: QualityControlOutput::from_files(&id, &files.qc)? ,
            ncbi_blast: match files.ncbi_blast { 
                Some(ref path) => String::new(), 
                None => {
                    log::error!("No NCBI BLAST file detected for: {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
        })
    }
}