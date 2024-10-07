use std::path::PathBuf;

use vircov::vircov::VircovSummary;

use crate::error::WorkflowError;
use crate::modules::pathogen::{KrakenReport, MetabuliReport};
use crate::utils::{
    get_file_by_name, get_file_component, FileComponent
};
use super::quality::{QualityControlFiles, QualityControlOutput};

#[derive(Debug, Clone)]
pub struct PathogenProfileFiles {
    kraken2: Option<PathBuf>,
    metabuli: Option<PathBuf>,
    sylph: Option<PathBuf>,
    kmcp: Option<PathBuf>,
    bracken: Option<PathBuf>,
    vircov: Option<PathBuf>
}
impl PathogenProfileFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {

        Ok(Self {
            vircov: get_file_by_name(&path, &id, ".alignment.tsv")?,
            kraken2: get_file_by_name(&path, &id, ".kraken2.report")?,
            bracken: get_file_by_name(&path, &id, ".bracken.report")?,
            metabuli: get_file_by_name(&path, &id, ".metabuli.report")?,
            sylph: get_file_by_name(&path, &id, ".sylph.mpa")?,
            kmcp: get_file_by_name(&path, &id, ".kmcp.profile")?,
        })

    }
}

#[derive(Debug, Clone)]
pub struct PathogenAssemblyFiles {

}
impl PathogenAssemblyFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self { })
    }
}


pub struct PathogenAssemblyOutput {
    
}
impl PathogenAssemblyOutput {
    pub fn from_files(id: &str, files: &PathogenProfileFiles) -> Result<Self, WorkflowError> {

        Ok(Self {
            
        })
    }
}

#[derive(Debug, Clone)]
pub struct PathogenFiles {
    pub qc: QualityControlFiles,
    pub profile: PathogenProfileFiles,
    pub assembly: PathogenAssemblyFiles
}
impl PathogenFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            qc: QualityControlFiles::from(path, id)?,
            profile: PathogenProfileFiles::from(path, id)?,
            assembly: PathogenAssemblyFiles::from(path, id)?
        })
    }
}

pub struct PathogenProfileOutput {
    pub id: String,
    pub vircov: Option<VircovSummary>,
    pub kraken2: Option<KrakenReport>,
    pub metabuli: Option<MetabuliReport>,

    // pub bracken: Option<BrackenProfile>,
    // pub sylph: Option<SylphProfile>,
    // pub kmcp: Option<KmcpProfile>
}
impl PathogenProfileOutput {
    pub fn from_files(id: &str, files: &PathogenProfileFiles) -> Result<Self, WorkflowError> {

        Ok(Self {
            id: id.to_string(),
            vircov: match files.vircov { 
                Some(ref path) => Some(VircovSummary::from_tsv(path, true)?), 
                None => None
            },
            kraken2: match files.kraken2 {
                Some(ref path) => Some(KrakenReport::from_report(path, true)?), 
                None => None
            },
            metabuli: match files.metabuli {
                Some(ref path) => Some(MetabuliReport::from_report(path, true)?), 
                None => None
            },
        })
    }
}

pub struct PathogenOutput {
    pub id: String,
    pub qc: QualityControlOutput,
    pub profile: PathogenProfileOutput,
    pub assembly: PathogenAssemblyOutput
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
            qc: QualityControlOutput::from_files(&id, &files.qc, background)?,
            profile: PathogenProfileOutput::from_files(&id, &files.profile)?,
            assembly: PathogenAssemblyOutput::from_files(&id, &files.profile)?
        })
    }
}

