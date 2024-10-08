use std::{fs::File, path::PathBuf, io::{Write, Read}};
use serde::{Serialize, Deserialize};
use crate::api::cerebro::model::ModelError;

use super::model::Pipeline;

type StageId = String;
type TowerId = String;

#[derive(Deserialize, Serialize, Debug)]
pub struct RegisterTowerSchema {
    pub id: TowerId,
    pub stage: StageId,
    pub date: String,
    pub name: String,
    pub location: String,
    pub last_ping: String,
    pub pipelines: Vec<Pipeline>,
}
impl RegisterTowerSchema {
    pub fn new(name: &str, location: &str, pipelines: Option<Vec<Pipeline>>) -> Self {
        
        let date = chrono::Utc::now().to_string();

        let pipelines = match pipelines {
            Some(pipelines) => pipelines,
            None => Vec::from([Pipeline::PanviralEnrichment, Pipeline::PathogenDetection, Pipeline::CultureIdentification])
        };

        Self {
            id: uuid::Uuid::new_v4().to_string(),
            date: date.clone(),
            name: name.to_string(),
            location: location.to_string(),
            last_ping: date,
            pipelines,
            stage: uuid::Uuid::new_v4().to_string()
        }
    }
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