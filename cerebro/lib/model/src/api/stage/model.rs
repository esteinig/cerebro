use std::{fs::File, path::PathBuf};
use std::io::{Read, Write};
use serde::{Serialize, Deserialize};

use crate::api::towers::model::Pipeline;
use crate::api::{cerebro::model::ModelError, files::model::SeaweedFile, towers::model::ProductionTower};


/*
========================
File system and storage
========================
*/

pub type FileId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StagedSample {
    pub id: String,
    pub date: String,
    pub run_id: String,
    pub sample_id: String,
    pub database: String,
    pub project: String,
    pub pipeline: Pipeline,
    pub tower: ProductionTower,
    pub files: Vec<SeaweedFile>
}
impl StagedSample {

    pub fn from_json(path: &PathBuf) -> Result<Self, ModelError> {
        let mut file = File::open(path)?;
        let mut json_data = String::new();
        file.read_to_string(&mut json_data)?;
        let schema = serde_json::from_str(&json_data).map_err(ModelError::JsonDeserialization)?;
        Ok(schema)
    }
    pub fn to_json(&self, path: &PathBuf) -> Result<(), ModelError> {
        let json_data = serde_json::to_string_pretty(&self).map_err(ModelError::JsonSerialization)?;
        let mut file = File::create(path)?;
        file.write_all(json_data.as_bytes())?;

        Ok(())
    }
}