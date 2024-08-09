use serde::{Serialize, Deserialize};
use super::schema::RegisterPipelineSchema;

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
    pub fn from_schema(schema: &RegisterPipelineSchema) -> Self {
        
        Self {
            id: schema.id.clone(),
            date: schema.date.clone(),
            name: schema.name.clone(),
            location: schema.location.clone(),
            last_ping: schema.last_ping.clone(),
            pipeline: schema.pipeline.clone(),
        }
    }
}

