use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}, path::{Path, PathBuf}};

use serde::{Deserialize, Serialize};
use taxonomy::ncbi;
use vircov::vircov::VircovRecord;

use crate::{error::WorkflowError, nextflow::{panviral::PanviralOutput, pathogen::PathogenDetectionOutput}, taxa::taxon::{Taxon, TaxonExtraction}, utils::read_tsv};


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Alignment {
    pub id: String,
    pub records: Vec<VircovRecord>
}
impl Alignment {
    pub fn from_panviral(
        output: &PanviralOutput,
    ) -> Result<Self, WorkflowError> {
       
        Ok(Self {
            id: output.id.clone(),
            records: output.vircov.records.clone()
        })

    }
    pub fn from_pathogen(
        output: &PathogenDetectionOutput,
    ) -> Result<Self, WorkflowError> {
       
        Ok(Self {
            id: output.id.clone(),
            records: output.profile.vircov.clone().map(|summary| summary.records).unwrap_or_default()
        })

    }
    pub fn from_tsv(&self, path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self {
            id: id.to_string(),
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
}

impl TaxonExtraction for Alignment {
    fn get_taxa(&self, taxonomy_directory: &PathBuf, strict: bool) -> Result<HashMap<String, Taxon>, WorkflowError> {
        let taxonomy = ncbi::load(taxonomy_directory)?;
    
        let mut taxa: HashMap<String, Taxon> = HashMap::new();
        for record in &self.records {
            let taxid = match &record.taxid {
                None => {
                    log::error!("Taxid not found for record with reference sequence: {}", record.reference);
                    return Err(WorkflowError::PanviralTaxidAnnotationMissing);
                },
                Some(taxid) => taxid,
            };
    
            if let Some(existing_taxon) = taxa.get_mut(taxid) {
                // If the taxid is already present, update its evidence
                existing_taxon.evidence.alignment.push(record.clone());
            } else {
                // Otherwise, create a new Taxon and insert it into the HashMap
                let mut taxon = match Taxon::from_taxid(taxid.clone(), &taxonomy, true) {
                    Err(err) => {
                        log::error!("Failed to find taxid '{}' in provided taxonomy", taxid);
                        if strict {
                            return Err(err);
                        } else {
                            log::warn!(
                                "Strict mode is not enabled - detection of taxid '{}' will be skipped",
                                taxid
                            );
                            continue;
                        }
                    },
                    Ok(taxon) => taxon,
                };
                taxon.evidence.alignment.push(record.clone());
                taxa.insert(taxid.clone(), taxon);
            }
        }
        Ok(taxa)
    }
}