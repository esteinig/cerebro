use serde::{Serialize, Deserialize};

use super::schema::RegisterPipelineSchema;

/*
========================
File system and storage
========================
*/


#[derive(Debug, Clone, Serialize, Deserialize, clap::ValueEnum)]
pub enum Pipeline {
    PathogenDetection,
    PanviralEnrichment,
    CultureIdentification
}
impl std::fmt::Display for Pipeline {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::PathogenDetection => write!(f, "PathogenDetection"),
            Self::PanviralEnrichment => write!(f, "PanviralEnrichment"),
            Self::CultureIdentification => write!(f, "CultureIdentification"),

        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProductionPipeline {
    pub id: String,
    pub date: String,
    pub name: String,
    pub location: String,
    pub last_ping: String,
    pub pipeline: Pipeline
}
impl ProductionPipeline {
    pub fn from_schema(register_production_pipeline_schema: &RegisterPipelineSchema) -> Self {
        Self {
            id: register_production_pipeline_schema.id.clone(),
            date: register_production_pipeline_schema.date.clone(),
            name: register_production_pipeline_schema.name.clone(),
            location: register_production_pipeline_schema.location.clone(),
            last_ping: register_production_pipeline_schema.last_ping.clone(),
            pipeline: register_production_pipeline_schema.pipeline.clone(),
            
        }
    }
}

