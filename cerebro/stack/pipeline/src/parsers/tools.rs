
use serde::{Serialize, Deserialize};
use std::io::{BufRead, BufReader};
use taxonomy::GeneralTaxonomy;
use std::collections::HashMap;
use std::path::PathBuf; 
use std::fs::File;

use crate::{
    taxon::{Taxon, TaxonomyWarning}, 
    record::{BlastLcaRecord, BlastRecord, DiamondRecord, Kraken2UniqRecord, VircovRecord, VircovScanRemapRecord}, 
    error::WorkflowError,
    virus::{AnnotationOptions, Annotation}
};

const WORKFLOW_RESULT_SPLIT: &str = "__";


fn get_meta_data_from_result_file(path: &PathBuf) -> Result<(String, String, String), WorkflowError> {

    let name = crate::utils::get_file_stem(path)?;
    let file_meta = name.split(WORKFLOW_RESULT_SPLIT).collect::<Vec<&str>>();
    let module = match file_meta.get(0) {
        Some(value) => value,
        None => return Err(WorkflowError::ResultFileIndexNotFound(name, 0.to_string()))
    }.to_string();
    let tool = match file_meta.get(1) {
        Some(value) => value,
        None => return Err(WorkflowError::ResultFileIndexNotFound(name, 1.to_string()))
    }.to_string();
    let db = match file_meta.get(2) {
        Some(value) => value,
        None => return Err(WorkflowError::ResultFileIndexNotFound(name, 2.to_string()))
    }.to_string();

    return Ok((module, tool, db))
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Blast {
    pub db: String,
    pub taxa: HashMap<String, Taxon>,
    pub records: Vec<BlastLcaRecord>,
    pub taxonomy_warnings: Vec<TaxonomyWarning>
}

impl Blast {
    ///
    /// 
    pub fn from(id: String, blast_report: &PathBuf, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {

        let (_, _, db) = get_meta_data_from_result_file(blast_report)?;
        let report_file = BufReader::new(File::open(blast_report)?);
                
        let mut records: Vec<BlastRecord> = vec![];
        for line in report_file.lines() {
            records.push(BlastRecord::from_str(line?, &id, &db)?)
        }

        let mut contigs: HashMap<String, Vec<BlastRecord>> = HashMap::new();
        for record in records {
            let contig = contigs.entry(record.qid.clone()).or_insert(vec![]);
            contig.push(record)
        }

        let mut lca_records: Vec<BlastLcaRecord> = vec![];
        let mut taxonomy_warnings: Vec<TaxonomyWarning> = vec![];
        for (_, blast_records) in contigs {
            
           let blast_lca_results = BlastLcaRecord::from_blast_record(&id, &db, blast_records, taxonomy);
            // At the moment we remove non-matching taxonomy queries from result - we need to completely synchronize databases with current taxonomy!
            // These are added to a rejection string record:
            if let Err(_) = blast_lca_results {
                taxonomy_warnings.push(TaxonomyWarning::new(&id, "Blast", &db, "BlastLcaRecord::from", format!("{}", blast_lca_results.unwrap_err())));
                continue;
            }
            lca_records.push(blast_lca_results?);

        }

        Ok(Self { db, taxa: HashMap::new(), records: lca_records, taxonomy_warnings })
    }
    pub fn get_taxa(&mut self, taxonomy: &GeneralTaxonomy) -> Result<(), WorkflowError>{
        if !self.taxa.is_empty() {
            return Err(WorkflowError::TaxAggregate("Blast".to_string()))
        }
        for record in &self.records {
            // Get the taxon for this record from its associated taxid
            let mut _taxon_result = Taxon::from_taxid(record.taxid.clone(), taxonomy, true);

            // At the moment we remove non-matching taxonomy queries from result - we need to completely synchronize databases with current taxonomy.
            if let Err(_) = _taxon_result {
                log::warn!("Could not find taxid `{}` in taxonomy! Entry is removed! Please ensure that databases and taxonomy are synchronized.", record.taxid);
                log::warn!("Record is: {} {} {}", record.id, record.taxid, record.as_row_field());
                self.taxonomy_warnings.push(TaxonomyWarning::new(&record.id, "Blast", &record.db, "Blast::get_taxa", format!("{} {}", record.taxid, record.as_row_field())));
                continue;
            }

            let mut _taxon = _taxon_result?;

            // If the taxon already exists in the collection, add the evidence to the existing taxon
            if let Some(taxon) = self.taxa.get_mut(&_taxon.taxid) {
                taxon.evidence.assembly.push(record.clone());
            } else {
                // ... otherwise add the evidence to the just created taxon and insert it into colletion
                _taxon.evidence.assembly.push(record.clone());
                self.taxa.insert(record.taxid.clone(), _taxon);
            }
        }
        Ok(())
    }
}


/// Blastn records and LCA processing
/// 
/// 
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Diamond {
    pub db: String,
    pub taxa: HashMap<String, Taxon>,
    pub records: Vec<BlastLcaRecord>,
    pub taxonomy_warnings: Vec<TaxonomyWarning>
}

impl Diamond {
    ///
    /// 
    pub fn from(id: String, diamond_report: &PathBuf, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {

        let (_, _, db) = get_meta_data_from_result_file(diamond_report)?;
        let report_file = BufReader::new(File::open(diamond_report)?);
                
        let mut records: Vec<DiamondRecord> = vec![];
        for line in report_file.lines() {
            records.push(DiamondRecord::from_str(line?, &id, &db)?)
        }

        let mut contigs: HashMap<String, Vec<DiamondRecord>> = HashMap::new();
        for record in records {
            let contig = contigs.entry(record.qid.clone()).or_insert(vec![]);
            contig.push(record)
        }

        let mut lca_records: Vec<BlastLcaRecord> = vec![];
        let mut taxonomy_warnings: Vec<TaxonomyWarning> = vec![];
        for (_, diamond_records) in contigs {
            
           let blast_lca_results = BlastLcaRecord::from_diamond_record(&id, &db, diamond_records, taxonomy);
            // At the moment we remove non-matching taxonomy queries from result - we need to completely synchronize databases with current taxonomy!
            // These are added to a warning record and added to the taxonomy warnings that are bubbled up and can be displayed in a table.
            if let Err(_) = blast_lca_results {
                taxonomy_warnings.push(TaxonomyWarning::new(&id, "Diamond", &db, "BlastLcaRecord::from", format!("{}", blast_lca_results.unwrap_err())));
                continue;
            }
            lca_records.push(blast_lca_results?);
        }

        Ok(Self { db, taxa: HashMap::new(), records: lca_records, taxonomy_warnings })
    }
    pub fn get_taxa(&mut self, taxonomy: &GeneralTaxonomy) -> Result<(), WorkflowError>{
        if !self.taxa.is_empty() {
            return Err(WorkflowError::TaxAggregate("Diamond".to_string()))
        }
        for record in &self.records {
            // Get the taxon for this record from its associated taxid
            let mut _taxon_result = Taxon::from_taxid(record.taxid.clone(), taxonomy, true);


            // At the moment we remove non-matching taxonomy queries from result - we need to completely synchronize databases with current taxonomy.
            if let Err(_) = _taxon_result {
                log::warn!("Could not find taxid `{}` in taxonomy! Entry is removed! Please ensure that databases and taxonomy are synchronized.", record.taxid);
                log::warn!("Record is: {} {} {}", record.id, record.taxid, record.as_row_field());
                self.taxonomy_warnings.push(TaxonomyWarning::new(&record.id, "Diamond", &record.db, "Diamond::get_taxa", format!("{} {}", record.taxid, record.as_row_field())));
                continue;
            }

            let mut _taxon = _taxon_result?;

            // If the taxon already exists in the collection, add the evidence to the existing taxon
            if let Some(taxon) = self.taxa.get_mut(&_taxon.taxid) {
                taxon.evidence.assembly.push(record.clone());
            } else {
                // ... otherwise add the evidence to the just created taxon and insert it into colletion
                _taxon.evidence.assembly.push(record.clone());
                self.taxa.insert(record.taxid.clone(), _taxon);
            }
        }
        Ok(())
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Kraken2Uniq {
    pub db: String,
    pub taxa: HashMap<String, Taxon>,
    pub records: Vec<Kraken2UniqRecord>,
    pub taxonomy_warnings: Vec<TaxonomyWarning>
}

impl Kraken2Uniq {
    pub fn from(id: String, kraken_report: &PathBuf) -> Result<Self, WorkflowError> {

        let (_, _, db) = get_meta_data_from_result_file(kraken_report)?;
        let report_file = BufReader::new(File::open(kraken_report)?);
        
        let mut records: Vec<Kraken2UniqRecord> = Vec::new();

        for line in report_file.lines() {
            let record: Kraken2UniqRecord = Kraken2UniqRecord::from_str(line?, id.clone(), db.clone())?;

            // Account for unclassified record with taxid="0" and taxname="Unclassified" which are not part of the NCBI Taxonomy
            // and only parse records with direct read assignments, taxonomy is adjusted with the get_taxa method. If not reading
            // direct assigned records, we would read the whole tree at all levels, but we adjust anyway so no need. It also means
            // we can aggregate at species level more easily further up the hierarchy.
            if !(record.taxid == "0" || record.taxname.to_lowercase() == "unclassified") && (record.reads_direct > 0) {
                records.push(record.clone());
            }
        }
        Ok(Self { db, taxa: HashMap::new(), records, taxonomy_warnings: Vec::new() })
    }
    pub fn get_taxa(&mut self, taxonomy: &GeneralTaxonomy) -> Result<(), WorkflowError>{
        if !self.taxa.is_empty() {
            return Err(WorkflowError::TaxAggregate("Kraken2Uniq".to_string()))
        }
        // Calling this twice will overwrite the current taxa with new taxons
        for record in &self.records {
            
            // Get the taxon for this record from its associated taxid
            let mut _taxon_result = Taxon::from_taxid(record.taxid.clone(), taxonomy, true);

            // At the moment we remove non-matching taxonomy queries from result - we need to completely synchronize databases with current taxonomy.
            if let Err(_) = _taxon_result {
                log::warn!("Could not find taxid `{}` in taxonomy! See the warning table output. Please ensure that databases and taxonomy are synchronized.", record.taxid);
                self.taxonomy_warnings.push(TaxonomyWarning::new(&record.id, "Kraken2Uniq", &record.db, "Kraken2Uniq::get_taxa", format!("{} {} {}", record.taxid, record.as_row_field(), record.taxname)));
                continue;
            }

            let mut _taxon = _taxon_result?;
            
            // If the taxon already exists in the collection, add the evidence to the existing taxon
            if let Some(taxon) = self.taxa.get_mut(&_taxon.taxid) {
                taxon.evidence.kmer.push(record.clone());
            } else {
                // ... otherwise add the evidence to the just created taxon and insert it into colletion
                _taxon.evidence.kmer.push(record.clone());
                self.taxa.insert(record.taxid.clone(), _taxon);
            }
        }
        Ok(())
    }

}

