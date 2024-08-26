use std::path::PathBuf;

use serde::{Serialize, Deserialize};

use crate::api::{cerebro::model::ModelError, pipelines::schema::RegisterPipelineSchema};

type PipelineId = String;

#[derive(Deserialize, Serialize, Debug)]
pub struct RegisterStagedSampleSchema {
    pub id: PipelineId,
    pub file_ids: Option<Vec<String>>,
    pub run_id: Option<String>,
}
impl RegisterStagedSampleSchema {
    pub fn new(id: &str, file_ids: Option<Vec<String>>, run_id: Option<String>) -> Self {
        Self {
            id: id.to_string(),
            file_ids,
            run_id,
        }
    }
    pub fn from_pipeline_json(path: &PathBuf, file_ids: Option<Vec<String>>, run_id: Option<String>) -> Result<Self, ModelError> {
        let schema = RegisterPipelineSchema::from_json(path)?;
        Ok(Self {
            id: schema.id,
            file_ids,
            run_id
        })

    }
}