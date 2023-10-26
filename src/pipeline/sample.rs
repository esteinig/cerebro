
use std::fs::File;
use anyhow::Result;
use std::path::PathBuf;
use itertools::Itertools;
use std::collections::HashMap;
use taxonomy::GeneralTaxonomy;
use std::io::{BufReader, Write};
use super::error::WorkflowError;
use super::output::WorkflowOutputs;
use serde::{Deserialize, Serialize};
use crate::pipeline::taxon::aggregate;
use crate::tools::modules::virus::AnnotationOptions;
use super::filters::TaxonFilterSettings;
use super::module::{QualityControlModule, KmerModule, AlignmentModule, AssemblyModule};
use super::taxon::{TaxonomyWarning, Taxon, TaxonEvidence};



// Read a filter settings object from JSON file
fn _read_filter_settings(file: &PathBuf) -> Result<TaxonFilterSettings, WorkflowError> {
    let reader = BufReader::new(File::open(file)?);
    let filter_settings: TaxonFilterSettings = serde_json::from_reader(reader).map_err(WorkflowError::JsonSerialization)?;
    Ok(filter_settings)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkflowSample {
    pub id: String,
    pub taxa: HashMap<String, Taxon>,
    pub qc_module: QualityControlModule,
    pub kmer_module: KmerModule,
    pub alignment_module: AlignmentModule,
    pub assembly_module: AssemblyModule,
    pub taxonomy_warnings: Vec<TaxonomyWarning>  // bubbled up taxonomy warnings
}
impl WorkflowSample {
    pub fn new(path: &PathBuf, id: &Option<String>, taxonomy: Option<PathBuf>, species_evidence: bool) -> Result<Self, WorkflowError> {

        log::info!("Parsing pipeline outputs in: {}", &path.display());
        let files = WorkflowOutputs::from_results(path)?;

        let id = match id {
            Some(value) => value.to_string(),
            None => crate::pipeline::utils::get_file_stem(&path)?
        };
        log::info!("ID assigned to sample: {}", &id);
        let qc_module = QualityControlModule::from(id.to_owned(), &files)?;

        match &taxonomy {
            Some(taxonomy) => {


                let taxonomy = taxonomy::ncbi::load(&taxonomy).map_err(|err| WorkflowError::TaxonomyNotParsed(err))?;

                // Kmer module instantiation
                let kmer_module = KmerModule::from(
                    id.to_owned(), files.kraken2uniq, &taxonomy
                )?;
        
                // Alignment module
                let alignment_module = AlignmentModule::from(
                    id.to_owned(), files.alignment, files.alignment_remap, &taxonomy, AnnotationOptions::virosaurus()
                )?;
        
                // Assembly module
                let assembly_module = AssemblyModule::from(
                    id.to_owned(), files.blastn, files.diamond, &taxonomy
                )?;
        
                let wfs = WorkflowSample::aggregate_modules(&mut Self { 
                    id: id.clone(), 
                    taxa: HashMap::new(), 
                    qc_module, 
                    kmer_module, 
                    alignment_module, 
                    assembly_module, 
                    taxonomy_warnings: Vec::new() 
                })?;
                
                if species_evidence {
                    log::info!("{:<10} - Aggregating taxa at species level", &id);
                    wfs.aggregate_species_evidence(&taxonomy)
                } else {
                    Ok(wfs)
                }
            },
            None => {
                
                log::warn!("No taxonomy provided for sample {} - processed only quality control data", id);

                Ok(Self {
                    id: id.clone(),
                    taxa: HashMap::new(),
                    qc_module: qc_module, 
                    kmer_module: KmerModule::new(), 
                    alignment_module: AlignmentModule::new(), 
                    assembly_module: AssemblyModule::new(), 
                    taxonomy_warnings: Vec::new()
                })
            }
        }
        

    }
    pub fn read_json(input: &PathBuf) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(input)?);
        serde_json::from_reader(reader).map_err(WorkflowError::JsonSerialization)
    }
    pub fn write_json(&self, file: &PathBuf) -> Result<(), WorkflowError> {
        let mut file = File::create(&file)?;
        let json_string = serde_json::to_string_pretty(&self).map_err(WorkflowError::JsonSerialization)?;
        write!(file, "{}", json_string)?;
        Ok(())
    }
    /// Aggregate taxa across modules
    fn aggregate_modules(&mut self) -> Result<Self, WorkflowError> {
        if !self.taxa.is_empty() {
            return Err(WorkflowError::TaxAggregate("WorfklowSample".to_string()))
        }
        let mut cc = self.clone();

        log::info!("Aggregating taxa from KmerModule");
        cc.taxa = aggregate(&mut cc.taxa, &self.kmer_module.taxa);

        log::info!("Aggregating taxa from AlignmentModule");
        cc.taxa = aggregate(&mut cc.taxa, &self.alignment_module.taxa);

        log::info!("Aggregating taxa from AssemblyModule");
        cc.taxa = aggregate(&mut cc.taxa, &self.assembly_module.taxa);

        cc.taxonomy_warnings.append(&mut self.kmer_module.taxonomy_warnings);
        cc.taxonomy_warnings.append(&mut self.alignment_module.taxonomy_warnings);
        cc.taxonomy_warnings.append(&mut self.assembly_module.taxonomy_warnings);

        // Aggregate and then update reads per million in evidence records for convenience
        cc = cc.update_rpm();
        cc = cc.update_bpm();

        Ok(cc)
    }
    /// Aggregate taxa and their evidence at this parent rank (upwards e.g. below species to species and none above)
    /// this will automatically remove all other taxa that do not have classifications at species level or below.
    pub fn aggregate_species_evidence(&self, taxonomy: &GeneralTaxonomy) -> Result<WorkflowSample, WorkflowError> {

        let mut parent_taxids: Vec<Option<String>> = Vec::new();
        for (_taxid, taxon) in &self.taxa {
            parent_taxids.push(taxon.level.species_taxid.clone())
        }

        let unique_taxids: Vec<String> = parent_taxids.iter().flatten().map(String::from).unique().collect();

        let mut new_taxa = HashMap::new();
        for taxid in unique_taxids {
            let taxid_evidence: &Vec<TaxonEvidence> = &self.taxa.iter().filter(|(_, taxon)| {
                taxon.level.species_taxid == Some(taxid.clone())
            }).map(|(_, taxon)| taxon.evidence.clone()).collect();

            let mut new_level_taxon = Taxon::from_taxid(taxid.clone(), taxonomy, true)?;

            new_level_taxon.evidence = TaxonEvidence::from(taxid_evidence);
            new_taxa.insert(taxid, new_level_taxon);

        }
        
        let mut cc = self.clone();
        cc.taxa = new_taxa;

        Ok(cc)

    }
    pub fn update_rpm(&self) -> Self {

        let mut cc = self.clone();
        for (_, taxon) in &mut cc.taxa {

            if let Some(fastp_data) = &self.qc_module.fastp {
                // Kmer evidence
                let rpm_kmer_records = taxon.evidence.kmer.clone().into_iter().map(|mut record| {
                    record.rpm = (record.reads as f64 / fastp_data.summary.before.reads as f64)*1_000_000 as f64;
                    record
                }).collect();
                taxon.evidence.kmer = rpm_kmer_records;

                // Alignment evidence
                let rpm_alignment_records = taxon.evidence.alignment.clone().into_iter().map(|mut record| {
                    record.scan_rpm = (record.scan_reads as f64 / fastp_data.summary.before.reads as f64)*1_000_000 as f64;  // using scan alignments as read counts! make sure this is ok
                    record
                }).collect();
                taxon.evidence.alignment = rpm_alignment_records;
            }
        }
        cc
    } 
    pub fn update_bpm(&self) -> Self {

        let mut cc = self.clone();
        for (_, taxon) in &mut cc.taxa {

            if let Some(fastp_data) = &self.qc_module.fastp {
                // Assembly evidence
                let bpm_assembly_records = taxon.evidence.assembly.clone().into_iter().map(|mut record| {
                    record.bpm = (record.length as f64 / fastp_data.summary.before.reads as f64)*1_000_000 as f64;
                    record
                }).collect();
                taxon.evidence.assembly = bpm_assembly_records;

            }
        }
        cc
    } 
}




