use serde::{Serialize, Deserialize};

use crate::api::{files::model::SeaweedFile, pipelines::model::ProductionPipeline};


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
    pub pipeline: ProductionPipeline,
    pub files: Vec<SeaweedFile>
}
