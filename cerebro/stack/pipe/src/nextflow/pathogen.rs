use std::collections::HashMap;
use std::path::PathBuf;

use taxonomy::ncbi;
use vircov::vircov::VircovSummary;

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::utils::{is_file_empty, read_tsv};
use crate::error::WorkflowError;
use crate::utils::{get_file_by_name, get_file_component, FileComponent};
use super::mag::{BlastTaxidMethod, MetagenomeAssemblyFiles, MetagenomeAssemblyOutput};
use super::quality::{QualityControlFiles, QualityControlOutput};

#[derive(Debug, Clone)]
pub struct PathogenProfileFiles {
    vircov: Option<PathBuf>,
    kraken2: Option<PathBuf>,
    metabuli: Option<PathBuf>,
    sylph: Option<PathBuf>,
    kmcp_reads: Option<PathBuf>,
    kmcp_abundance: Option<PathBuf>,
    ganon_reads: Option<PathBuf>,
    ganon_abundance: Option<PathBuf>,
    bracken: Option<PathBuf>,
}
impl PathogenProfileFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {

        Ok(Self {
            vircov: get_file_by_name(&path, &id, ".vircov.tsv")?,
            kraken2: get_file_by_name(&path, &id, ".kraken2.reads.report")?,
            metabuli: get_file_by_name(&path, &id, ".metabuli.reads.report")?,
            ganon_reads: get_file_by_name(&path, &id, ".ganon.reads.report")?,
            bracken: get_file_by_name(&path, &id, ".bracken.abundance.report")?,
            kmcp_reads: get_file_by_name(&path, &id, ".kmcp.reads.report")?,
            kmcp_abundance: get_file_by_name(&path, &id, ".kmcp.abundance.report")?,
            ganon_abundance: get_file_by_name(&path, &id, ".ganon.abundance.report")?,
            sylph: get_file_by_name(&path, &id, ".sylph.abundance.report")?,
        })

    }
}

#[derive(Debug, Clone)]
pub struct PathogenFiles {
    pub qc: QualityControlFiles,
    pub profile: PathogenProfileFiles,
    pub mag: MetagenomeAssemblyFiles
}
impl PathogenFiles {
    pub fn from(path: &PathBuf, path_qc: Option<PathBuf>, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            qc: QualityControlFiles::from(match path_qc { Some(ref p) => p, None => path}, id)?,
            profile: PathogenProfileFiles::from(path, id)?,
            mag: MetagenomeAssemblyFiles::from(&path, &id)?
        })
    }
}

pub struct PathogenOutput {
    pub id: String,
    pub vircov: Option<VircovSummary>,
    pub kraken2: Option<KrakenReport>,
    pub metabuli: Option<MetabuliReport>,
    pub bracken: Option<BrackenReport>,
    pub kmcp_reads: Option<KmcpReadsReport>,
    pub kmcp_abundance: Option<KmcpAbundanceReport>,
    pub sylph: Option<SylphReport>,
    pub ganon_reads: Option<GanonReadsReport>,
    pub ganon_abundance: Option<GanonAbundanceReport>,
}
impl PathogenOutput {
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
            kmcp_reads: match files.kmcp_reads {
                Some(ref path) => Some(KmcpReadsReport::from_report(path, &id)?), 
                None => None
            },
            kmcp_abundance: match files.kmcp_abundance {
                Some(ref path) => Some(KmcpAbundanceReport::from_report(path, &id)?), 
                None => None
            },
            sylph: match files.sylph {
                Some(ref path) => Some(SylphReport::from_report(path, &id)?), 
                None => None
            },
            ganon_reads: match files.ganon_reads {
                Some(ref path) => Some(GanonReadsReport::from_report(path, &id)?), 
                None => None
            },
            ganon_abundance: match files.ganon_abundance {
                Some(ref path) => Some(GanonAbundanceReport::from_report(path, &id)?), 
                None => None
            },
        })
    }
}

pub struct PathogenDetectionOutput {
    pub id: String,
    pub qc: QualityControlOutput,
    pub profile: PathogenOutput,
    pub assembly: MetagenomeAssemblyOutput
}
impl PathogenDetectionOutput {

    pub fn from(path: &PathBuf, path_qc: Option<PathBuf>, id: Option<String>, taxonomy: Option<PathBuf>, blast_taxid: BlastTaxidMethod) -> Result<Self, WorkflowError> {

        let taxonomy = match taxonomy {
            Some(path) => Some(ncbi::load(&path)?),
            None => {
                match blast_taxid {
                    BlastTaxidMethod::LCA => return Err(WorkflowError::TaxonomyNotProvided),
                    BlastTaxidMethod::HighestBitscore => None
                }
                
            }
        };

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let files = PathogenFiles::from(&path, path_qc, &id)?;
        
        Ok(Self{
            id: id.to_string(),
            qc: QualityControlOutput::from_files(&id, &files.qc)?,
            profile: PathogenOutput::from_files(&id, &files.profile)?,
            assembly: MetagenomeAssemblyOutput::from_files(&id, &files.mag, taxonomy.as_ref(), blast_taxid)?
        })
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
pub struct KmcpAbundanceReportRecord {
    pub taxid: String,
    pub tax_level: String,
    pub taxid_lineage: String,
    pub taxname_lineage: String,
    pub abundance: f64,
}


#[derive(Serialize, Deserialize, Clone)]
pub struct KmcpReadsReportRecord {
    pub taxid: String,
    pub rank: String,
    pub taxname: String,
    pub reads: u64,
}


#[derive(Serialize, Deserialize)]
pub struct GanonReportRecord {
    pub tax_level: String,
    pub taxid: String,
    pub taxid_lineage: String,
    pub taxname: String,
    pub unique: u64,
    pub shared: u64,
    pub children: u64,
    pub cumulative: u64,
    pub cumulative_percent: f64
}

#[derive(Serialize, Deserialize)]
pub struct SylphReportRecord {
    pub clade_name: String,
    pub relative_abundance: f64,
    pub sequence_abundance: f64
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
            records: if is_file_empty(&path)? { Vec::new() } else { read_tsv(path, false, false)? }     
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
            records: if is_file_empty(&path)? { Vec::new() } else { read_tsv(path, false, false)? }    
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
            records: if is_file_empty(&path)? { Vec::new() } else { read_tsv(path, false, true)? }    
        })
    }
}

#[derive(Serialize, Deserialize)]
pub struct KmcpAbundanceReport {
    pub id: String,
    pub path: PathBuf,
    pub records: Vec<KmcpAbundanceReportRecord>
}
impl KmcpAbundanceReport {
    pub fn from_report(path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self { 
            id: id.to_string(), 
            path: path.to_path_buf(), 
            records: if is_file_empty(&path)? { Vec::new() } else { read_tsv(path, false, false)? }   
        })
    }
}

#[derive(Serialize, Deserialize)]
pub struct KmcpReadsReport {
    pub id: String,
    pub path: PathBuf,
    pub records: Vec<KmcpReadsReportRecord>
}
impl KmcpReadsReport {
    pub fn from_report(path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self { 
            id: id.to_string(), 
            path: path.to_path_buf(), 
            records: if is_file_empty(&path)? { Vec::new() } else { read_tsv(path, false, true)? }  
        })
    }
    
    pub fn get_taxid_report(&self) -> Result<KmcpReadsReport, WorkflowError> {
        let mut grouped_records: HashMap<String, KmcpReadsReportRecord> = HashMap::new();

        for record in &self.records {
            let taxid = &record.taxid;
            if let Some(existing_record) = grouped_records.get_mut(taxid) {
                // Sum the reads
                existing_record.reads += record.reads;
            } else {
                // Insert a clone of the record
                grouped_records.insert(taxid.clone(), record.clone());
            }
        }

        // Collect the grouped records into a vector
        let new_records: Vec<KmcpReadsReportRecord> = grouped_records.into_values().collect();

        Ok(KmcpReadsReport {
            id: self.id.clone(),
            path: self.path.clone(),
            records: new_records,
        })
    }
}


#[derive(Serialize, Deserialize)]
pub struct GanonReadsReport {
    pub id: String,
    pub path: PathBuf,
    pub records: Vec<GanonReportRecord>
}
impl GanonReadsReport {
    pub fn from_report(path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self { 
            id: id.to_string(), 
            path: path.to_path_buf(), 
            records: if is_file_empty(&path)? { Vec::new() } else { read_tsv(path, false, false)?  }  
        })
    }
}


#[derive(Serialize, Deserialize)]
pub struct GanonAbundanceReport {
    pub id: String,
    pub path: PathBuf,
    pub records: Vec<GanonReportRecord>
}
impl GanonAbundanceReport {
    pub fn from_report(path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self { 
            id: id.to_string(), 
            path: path.to_path_buf(), 
            records: if is_file_empty(&path)? { Vec::new() } else {  read_tsv(path, false, false)? }  
        })
    }
}


#[derive(Serialize, Deserialize)]
pub struct SylphReport {
    pub id: String,
    pub path: PathBuf,
    pub records: Vec<SylphReportRecord>
}
impl SylphReport {
    pub fn from_report(path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self { 
            id: id.to_string(), 
            path: path.to_path_buf(), 
            records: if is_file_empty(&path)? { Vec::new() } else { read_tsv(path, false, true)? } 
        })
    }
}