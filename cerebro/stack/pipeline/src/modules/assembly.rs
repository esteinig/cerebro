use std::{fs::File, io::{BufReader, BufWriter}, path::Path};

use serde::{Deserialize, Serialize};

use crate::{error::WorkflowError, nextflow::{mag::MetagenomeAssemblyOutput, pathogen::PathogenDetectionOutput}, utils::read_tsv};


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ContigRecord {
    pub id: String
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MetagenomeAssembly {
    pub id: String,
    pub records: Vec<ContigRecord>
}
impl MetagenomeAssembly {
    pub fn from_mag(
        output: &MetagenomeAssemblyOutput,
    ) -> Result<Self, WorkflowError> {
        
        Ok(Self {
            id: output.id.clone(),
            records: Vec::new()
        })

    }
    pub fn from_pathogen(
        output: &PathogenDetectionOutput,
    ) -> Result<Self, WorkflowError> {
        
        Ok(Self {
            id: output.id.clone(),
            records: Vec::new()
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
