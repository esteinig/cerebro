use std::path::PathBuf;

use vircov::vircov::VircovSummary;

use core::fmt;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::utils::read_tsv;
use crate::error::WorkflowError;
use crate::utils::{get_file_by_name, get_file_component, FileComponent};
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
            kmcp: get_file_by_name(&path, &id, ".kmcp.report")?,
            sylph: get_file_by_name(&path, &id, ".sylph.report")?,
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
    pub bracken: Option<BrackenReport>,
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
                Some(ref path) => Some(KrakenReport::from_report(path, &id)?), 
                None => None
            },
            bracken: match files.bracken {
                Some(ref path) => Some(BrackenReport::from_report(path, &id)?), 
                None => None
            },
            metabuli: match files.metabuli {
                Some(ref path) => Some(MetabuliReport::from_report(path, &id)?), 
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



/// Enum representing the available classifiers.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Classifier {
    #[serde(rename="kraken2")]
    Kraken2,
    #[serde(rename="metabuli")]
    Metabuli,
}


/// Enum representing the available profilers.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Profiler {
    #[serde(rename="sylph")]
    Sylph,
    #[serde(rename="kmcp")]
    Kmcp,
    #[serde(rename="bracken")]
    Bracken,

}

/// Enum representing the available aligners.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum Aligner {
    #[serde(rename="bowtie2")]
    Bowtie2,
    #[serde(rename="minimap2")]
    Minimap2,
    #[serde(rename="strobealign")]
    Strobealign,
}


impl Aligner {
    // Used for identification of pre-built-indices
    pub fn short_name(&self) -> &str {
        match self {
            Aligner::Bowtie2 => "bt2",
            Aligner::Minimap2 => "mm2",
            Aligner::Strobealign => "sti",
        }
    }
}
impl fmt::Display for Aligner {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Aligner::Bowtie2 => write!(f, "bowtie2"),
            Aligner::Minimap2 => write!(f, "minimap2"),
            Aligner::Strobealign => write!(f, "strobealign"),
        }
    }
}


impl Classifier {
    // Used for identification of pre-built-indices
    pub fn short_name(&self) -> &str {
        match self {
            Classifier::Kraken2 => "kraken2",
            Classifier::Metabuli => "metabuli",
        }
    }
}
impl fmt::Display for Classifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Classifier::Kraken2 => write!(f, "kraken2"),
            Classifier::Metabuli => write!(f, "metabuli"),
        }
    }
}


impl Profiler {
    // Used for identification of pre-built-indices
    pub fn short_name(&self) -> &str {
        match self {
            Profiler::Sylph => "sylph",
            Profiler::Kmcp => "kmcp",
            Profiler::Bracken => "bracken",
        }
    }
}
impl fmt::Display for Profiler {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Profiler::Sylph => write!(f, "sylph"),
            Profiler::Bracken => write!(f, "bracken"),
            Profiler::Kmcp  =>  write!(f, "kmcp")
        }
    }
}


#[derive(Serialize, Deserialize)]
pub struct KrakenReportRecord {
    pub percent: String,
    pub reads: u64,
    pub reads_direct: u64,
    pub tax_level: String,
    pub taxid: String,
    pub taxname: String
}

#[derive(Serialize, Deserialize)]
pub struct MetabuliReportRecord {
    pub percent: String,
    pub reads: u64,
    pub reads_direct: u64,
    pub tax_level: String,
    pub taxid: String,
    pub taxname: String
}

#[derive(Serialize, Deserialize)]
pub struct BrackenReportRecord {
    #[serde(rename(deserialize = "name"))]
    pub taxname: String,
    #[serde(rename(deserialize = "taxonomy_id"))]
    pub taxid: String,
    #[serde(rename(deserialize = "taxonomy_lvl"))]
    pub tax_level: String,
    #[serde(rename(deserialize = "kraken_assigned_reads"))]
    pub kraken_reads: u64,
    #[serde(rename(deserialize = "added_reads"))]
    pub bracken_reads: u64,
    #[serde(rename(deserialize = "new_est_reads"))]
    pub reads: u64,
    #[serde(rename(deserialize = "fraction_total_reads"))]
    pub fraction: f64
}


#[derive(Serialize, Deserialize)]
pub struct KrakenReport {
    pub id: String,
    pub path: PathBuf,
    pub records: Vec<KrakenReportRecord>
}
impl KrakenReport {
    pub fn from_report(path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self { 
            id: id.to_string(), 
            path: path.to_path_buf(), 
            records: read_tsv(path, false, false)? 
        })
    }
}

#[derive(Serialize, Deserialize)]
pub struct MetabuliReport {
    pub id: String,
    pub path: PathBuf,
    pub records: Vec<MetabuliReportRecord>
}
impl MetabuliReport {
    pub fn from_report(path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self { 
            id: id.to_string(), 
            path: path.to_path_buf(), 
            records: read_tsv(path, false, false)? 
        })
    }
}

#[derive(Serialize, Deserialize)]
pub struct BrackenReport {
    pub id: String,
    pub path: PathBuf,
    pub records: Vec<BrackenReportRecord>
}
impl BrackenReport {
    pub fn from_report(path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self { 
            id: id.to_string(), 
            path: path.to_path_buf(), 
            records: read_tsv(path, false, true)? 
        })
    }
}

