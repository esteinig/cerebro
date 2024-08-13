use serde::{Serialize, Deserialize};
use crate::api::pipelines::model::Pipeline;

#[derive(Deserialize, Serialize, Debug)]
pub struct RegisterStagedSampleSchema {
    pub run_id: String,
    pub database: String,
    pub project: String,
    pub pipeline: Pipeline
}
