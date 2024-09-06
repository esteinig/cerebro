use std::path::PathBuf;

use crate::error::WorkflowError;
use crate::utils::{
    get_file_component, 
    FileComponent
};
use super::quality::{QualityControlFiles, QualityControlOutput};

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

pub struct PathogenOutput {
    pub id: String,
    pub qc: QualityControlOutput,
    
}
impl PathogenOutput {

    pub fn from(path: &PathBuf, id: Option<String>, background: bool) -> Result<Self, WorkflowError> {

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let files = PathogenFiles::from(&path, &id)?;
        
        Ok(Self{
            id: id.to_string(),
            qc: QualityControlOutput::from_files(&id, &files.qc, background)? 
        })
    }
}

