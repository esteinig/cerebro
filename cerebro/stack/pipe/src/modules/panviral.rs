use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}, path::{Path, PathBuf}};

use serde::{Deserialize, Serialize};
use taxonomy::ncbi;
use vircov::vircov::VircovRecord;

use crate::{error::WorkflowError, nextflow::panviral::PanviralOutput, taxa::taxon::Taxon, utils::read_tsv};

use super::quality::QualityControl;


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Panviral {
    pub id: String,
    pub paired_end: bool,
    pub records: Vec<VircovRecord>
}
impl Panviral {
    pub fn from_panviral(
        output: &PanviralOutput,
        quality: &QualityControl,
        paired_end: bool,
    ) -> Result<Self, WorkflowError> {
       
        Ok(Self {
            id: output.id.clone(),
            paired_end,
            records: output.vircov.records.clone()
        })

    }
    pub fn from_tsv(&self, path: &Path, id: &str, paired_end: bool) -> Result<Self, WorkflowError> {
        Ok(Self {
            id: id.to_string(),
            paired_end,
            records: read_tsv(path, false, true)?
        })
    }
    pub fn to_json(&self, path: &Path) -> Result<(), WorkflowError> {
        let writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }
    pub fn from_json(path: &Path) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(path)?);
        let pathogen = serde_json::from_reader(reader)?;
        Ok(pathogen)
    }
    // Uses the provided taxonomy to create the taxon structs for the database model
    pub fn get_taxa(&self, taxonomy_directory: &PathBuf, strict: bool) -> Result<HashMap<String, Taxon>, WorkflowError> {

        let taxonomy = ncbi::load(taxonomy_directory)?;

        let mut taxa = HashMap::new();
        for record in &self.records {

            let taxid = match &record.taxid {
                None => {
                    log::error!("Taxid not found for record with reference sequence: {}", record.reference);
                    return Err(WorkflowError::PanviralTaxidAnnotationMissing)
                },
                Some(taxid) => taxid
            };

            let mut taxon = match Taxon::from_taxid(taxid.clone(), &taxonomy, true) {
                Err(err) => {
                    log::error!("Failed to find taxid '{}' in provided taxonomy", taxid);
                    if strict {
                        return Err(err)
                    } else {
                        log::warn!("Strict mode is not enabled - detection of taxid '{}' will be skipped", taxid);
                        continue;
                    }
                },
                Ok(taxon) => taxon,
            };
            taxon.evidence.alignment.push(record.clone());
            taxa.insert(taxid.clone(), taxon);
        }
        Ok(taxa)
    }
}