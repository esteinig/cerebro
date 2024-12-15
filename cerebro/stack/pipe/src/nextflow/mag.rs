use std::collections::HashMap;
use std::fs::File;
use std::{collections::HashSet, path::Path};
use std::path::PathBuf;
use serde::{Deserialize, Serialize};
use taxonomy::{GeneralTaxonomy, Taxonomy};
use itertools::Itertools;
use crate::utils::{
    get_file_by_name, get_file_component, read_tsv, FileComponent
};

use std::io::{BufRead, BufReader};
use crate::error::WorkflowError;



#[derive(Debug, Clone)]
pub struct MagFiles {
    pub ncbi_blast: Option<PathBuf>,
}
impl MagFiles {
    pub fn from(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            ncbi_blast: get_file_by_name(&path, &id, ".blast.assembly.tsv")?,
        })
    }
}


#[derive(Debug, Clone)]
pub struct MagOutput {
    pub id: String,
    pub ncbi_blast: Vec<ContigBlastRecord>,
    
}
impl MagOutput {

    pub fn from(id: Option<String>, path: &PathBuf, taxonomy: Option<&GeneralTaxonomy>, blast_taxid: BlastTaxidMethod) -> Result<Self, WorkflowError> {

        let id = match id {
            Some(id) => id,
            None => get_file_component(&path, FileComponent::FileName)?
        };

        let files = MagFiles::from(&path, &id)?; 
        
        Ok(Self{
            id: id.to_string(),
            ncbi_blast: match files.ncbi_blast { 
                Some(ref path) => Self::parse_blast_records(&id, path, taxonomy, blast_taxid)?, 
                None => {
                    log::error!("No NCBI BLAST file detected for: {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
        })
    }
    pub fn from_files(id: &str, files: &MagFiles, taxonomy: Option<&GeneralTaxonomy>, blast_taxid: BlastTaxidMethod) -> Result<Self, WorkflowError> {
        
        Ok(Self{
            id: id.to_string(),
            ncbi_blast: match files.ncbi_blast { 
                Some(ref path) => Self::parse_blast_records(&id, path, taxonomy, blast_taxid)?, 
                None => {
                    log::error!("No NCBI BLAST file detected for: {id}");
                    return Err(WorkflowError::PipelineOutputNotFound)
                }
            },
        })
    }
    pub fn parse_blast_records(id: &str, path: &PathBuf, taxonomy: Option<&GeneralTaxonomy>, blast_taxid: BlastTaxidMethod) -> Result<Vec<ContigBlastRecord>, WorkflowError> {
    
        let records: Vec<BlastRecord> = read_tsv(&path, false, false)?;

        // Group records by qseqid
        let mut grouped_records: HashMap<String, Vec<BlastRecord>> = HashMap::new();
        for record in records {
            grouped_records
                .entry(record.qid.clone())
                .or_insert_with(Vec::new)
                .push(record);
        }


        // Construct ContigBlastRecord for each group
        let mut contig_records = Vec::new();
        for (_, records) in grouped_records {
            let contig_record = ContigBlastRecord::from_blast_records(id, &records, taxonomy, blast_taxid.clone())?;
            contig_records.push(contig_record);
        }

        Ok(contig_records)
    }
}



// qseqid qlen qstart qend sseqid slen sstart send length nident pident evalue bitscore staxid ssciname stitle
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlastRecord {
    pub qid: String,
    pub qlen: u64,
    pub qstart: u64,
    pub qend: u64,
    pub sid: String,
    pub slen: u64,
    pub sstart: u64,
    pub send: u64,
    pub length: u64,
    pub nident: u64,
    pub pident: f64,
    pub evalue: f64,
    pub bitscore: f64,
    pub taxid: String,
    pub taxname: String,
    pub title: String
}

/// Enumeration for how the final taxonomic identifier is obtained
#[derive(Debug, Clone)]
pub enum BlastTaxidMethod {
    LCA,             // Use Lowest Common Ancestor approach for taxid
    HighestBitscore, // Use the taxid of the record with the highest bitscore
}


// qseqid qlen qstart qend sseqid slen sstart send length nident pident evalue bitscore staxids sscinames stitle 
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct DiamondRecord {
    pub id: String,
    pub db: String,
    pub qid: String,
    pub qlen: u64,
    pub qstart: u64,
    pub qend: u64,
    pub sid: String,
    pub slen: u64,
    pub sstart: u64,
    pub send: u64,
    pub length: u64,
    pub nident: u64,
    pub pident: f64,
    pub evalue: f64,
    pub bitscore: f64,
    pub taxids: Vec<String>,
    pub taxnames: Vec<String>,
    pub title: String
}
 

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ContigBlastRecord {
    pub id: String,
    pub length: u64,    
    pub alignment: u64,
    pub coverage: f64,         
    pub identity: f64,
    pub evalue: f64,
    pub bitscore: f64,
    pub taxid: String,
    pub title: String,
    pub reference: String,
    pub reference_length: u64,
    pub read_coverage: f64,           // if present in contig name e.g. from spades output
    pub bpm: f64,
}
impl ContigBlastRecord {
    pub fn from_blast_records(id: &str, records: &Vec<BlastRecord>, taxonomy: Option<&GeneralTaxonomy>, method: BlastTaxidMethod) -> Result<Self, WorkflowError> {
        
        // Get contig query length and name
        let (record_length, record_name) = match records.get(0) {
            Some(record) => (record.qlen, record.qid.clone()),
            None => return Err(WorkflowError::BlastLcaRecordExtraction)
        };

        // Extraction of coverage from MetaSpades assembly output
        let record_name_components: Vec<&str> = record_name.split("_").collect();
        let read_coverage = match record_name_components.contains(&"cov") {
            true => match record_name_components.last(){
                Some(cov) => cov.parse::<f64>().map_err(WorkflowError::BlastLcaFloatFieldConversion)?,
                None => return Err(WorkflowError::BlastLcaCovExtraction)
            },
            false => 0.0
        };
        

        let (record, taxid) = match method {
            BlastTaxidMethod::HighestBitscore => {
                match records.into_iter().max_by(|a, b| a.bitscore.total_cmp(&b.bitscore)){
                    Some(record) => (record.clone(), record.taxid.to_owned()),
                    None => return Err(WorkflowError::BlastLcaRecordExtractionEvalue)
                }
            },
            BlastTaxidMethod::LCA => {

                let taxonomy = match taxonomy {
                    Some(taxonomy) => taxonomy,
                    None => return Err(WorkflowError::TaxonomyNotProvided)
                };

                // Find the LCA of all BlastRecords for this contig
                let taxids_unique: Vec<String> = records.iter().map(|x| x.taxid.to_string()).unique().collect();

                let (failed_taxids, lca_taxid) = get_lca_taxid(taxids_unique.clone(), taxonomy)?;
                
                let mut no_lca_taxid = false;
                if !failed_taxids.is_empty() {
                    let proportion_failed = failed_taxids.len() as f64 / taxids_unique.len() as f64;
                    if proportion_failed >= 0.5 {
                        log::warn!("Failed to detect a significant proportion of taxids in the reference taxonomy from all hits for contig {}: {:?} ({:.2} %). LCA classification for this contig is ignored.", &record_name, &failed_taxids, proportion_failed*100.0);
                        no_lca_taxid = true;
                    }
                };

                match records.into_iter().max_by(|a, b| a.bitscore.total_cmp(&b.bitscore)){
                    Some(record) => (record.clone(), if no_lca_taxid { record.taxid.to_owned() } else { lca_taxid } ),
                    None => return Err(WorkflowError::BlastLcaRecordExtractionEvalue)
                }
            }
        };
        
        Ok(Self { 
            id: id.to_string(),
            length: record_length,
            alignment: record.length,
            coverage: (record.length as f64/record.qlen as f64)*100.0,
            identity: record.pident,
            evalue: record.evalue,
            bitscore: record.bitscore,
            taxid: taxid.to_string(),
            read_coverage,
            reference: record.sid,
            reference_length: record.slen,
            title: record.title,
            bpm: 0.0
        })
    }
    pub fn from_diamond_records(id: &str, records: Vec<DiamondRecord>, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {
        let (record_length, record_name) = match records.get(0) {
            Some(record) => (record.qlen, record.qid.clone()),
            None => return Err(WorkflowError::BlastLcaRecordExtraction)
        };

        // Extraction of coverage from MetaSpades assembly output
        let record_name_components: Vec<&str> = record_name.split("_").collect();
        let read_coverage = match record_name_components.contains(&"cov") {
            true => match record_name_components.last(){
                Some(cov) => cov.parse::<f64>().map_err(WorkflowError::BlastLcaFloatFieldConversion)?,
                None => return Err(WorkflowError::BlastLcaCovExtraction)
            },
            false => 0.0
        };
        
        let mut taxids: Vec<String> = Vec::new();
        for record in &records {
            for taxid in &record.taxids {
                taxids.push(taxid.to_string())
            }
        };
        let taxids_unique: Vec<String> = taxids.iter().unique().map(|x| x.to_string()).collect();

        // Find the LCA of all hits against this contig
        let (failed_taxids, lca_taxid) = get_lca_taxid(taxids_unique.clone(), taxonomy)?;
       
        if !failed_taxids.is_empty() {
            let proportion_failed = failed_taxids.len() as f64 / taxids_unique.len() as f64;
            if proportion_failed >= 0.5 {
                log::warn!("Failed to detect a significant proportion (>= 50%) of unique taxids from all hits for contig {}: {:?} ({:.2} %) in the reference taxonomy. LCA classification for this contig is ignored.", &record_name, &failed_taxids, proportion_failed*100.0);
                return Err(WorkflowError::BlastLcaTaxidsNotFound)
            }
        }

        let highest_bitscore = match records.into_iter().max_by(|a, b| a.bitscore.total_cmp(&b.bitscore)){
            Some(record) => record,
            None => return Err(WorkflowError::BlastLcaRecordExtractionEvalue)
        };
        
        Ok(Self { 
            id: id.to_string(),
            length: record_length,
            alignment: highest_bitscore.length,
            coverage: (highest_bitscore.length as f64/highest_bitscore.qlen as f64)*100.0,
            identity: highest_bitscore.pident,
            evalue: highest_bitscore.evalue,
            bitscore: highest_bitscore.bitscore,
            taxid: lca_taxid.to_string(),
            read_coverage,
            reference: highest_bitscore.sid,
            reference_length: highest_bitscore.slen,
            title: highest_bitscore.title,
            bpm: 0.
        } )
    }
    
}



// Compute the least common ancestor from a vector of unique taxids
pub fn get_lca_taxid(taxids_unique: Vec<String>, taxonomy: &GeneralTaxonomy) -> Result<(Vec<String>, String), WorkflowError> {

    let mut failed_taxids = Vec::new();

    match taxids_unique.len() {
        0 => return Err(WorkflowError::BlastLcaRecordExtraction),
        1 => {
            // If there is just one taxid, simply return it
            let taxid_str = taxids_unique[0].as_str();
            // Check that the taxid exists in the taxonomy, this is a safety check for later taxonomy operations 
            match taxonomy.rank(taxid_str) {
                Err(_) => {
                    failed_taxids.push(taxid_str.to_owned());
                    return Ok((failed_taxids, taxid_str.to_owned()))
                },
                Ok(_) => return Ok((failed_taxids, taxid_str.to_owned()))
            }
        },
        _ => {


            // LCA by getting lineages in ordered form
            let mut lineages: Vec<Vec<&str>> = vec![];
            for taxid in taxids_unique.iter() {
                let taxid_str = taxid.as_str();
                match taxonomy.lineage(taxid_str) {
                    // If we can't get the lineage for this taxid, 
                    // add it to the list of failed taxids and continue
                    Err(_) => failed_taxids.push(taxid_str.to_owned()),
                    Ok(lineage) => lineages.push(lineage)
                }
            };
            // Transform lineages into sets
            let lineage_sets: Vec<HashSet<&str>> = lineages.clone().into_iter().map(|x| HashSet::from_iter(x)).collect();
            // Take the set intersection of lineages
            let common_lineages: HashSet<&str> = if let Some((first, rest)) = lineage_sets.split_first() {
                rest.iter().fold(first.clone(), |acc, i| {
                    acc.intersection(i).copied().collect()
                })
            } else {
                HashSet::new()
            };
            // For each taxid in order from the first lineage  check if its contained in the set intersection 
            // of lineages and return it, otherwise assign root ('root')
            let mut lca_taxid = "root";
            if let Some(first) = lineages.first() {
                for taxid in first {
                    if common_lineages.contains(taxid) {
                        lca_taxid = taxid;
                        break;
                    }
                }
            }
            Ok((failed_taxids, lca_taxid.to_string()))
        }
    }
}
