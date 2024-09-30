use std::path::PathBuf;

use serde::{Serialize, Deserialize};

use crate::api::{cerebro::model::ModelError, towers::{model::Pipeline, schema::RegisterTowerSchema}};

type TowerId = String;

#[derive(Deserialize, Serialize, Debug)]
pub struct RegisterStagedSampleSchema {
    pub tower_id: TowerId,
    pub pipeline: Pipeline,
    pub file_ids: Option<Vec<String>>,
    pub run_id: Option<String>,
}
impl RegisterStagedSampleSchema {
    pub fn new(tower_id: &str, pipeline: Pipeline, file_ids: Option<Vec<String>>, run_id: Option<String>) -> Self {
        Self {
            tower_id: tower_id.to_string(),
            pipeline,
            file_ids,
            run_id,
        }
    }
    pub fn from_tower_json(path: &PathBuf, pipeline: Pipeline, file_ids: Option<Vec<String>>, run_id: Option<String>) -> Result<Self, ModelError> {
        let schema = RegisterTowerSchema::from_json(path)?;
        Ok(Self {
            tower_id: schema.id,
            pipeline,
            file_ids,
            run_id
        })

    }
}