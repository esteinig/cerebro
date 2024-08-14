use serde::{Serialize, Deserialize};
use super::schema::RegisterPipelineSchema;

#[derive(Debug, Clone, Serialize, Deserialize, clap::ValueEnum)]
pub enum Pipeline {
    #[serde(rename="pathogen-detection")]
    PathogenDetection,
    #[serde(rename="panviral-enrichment")]
    PanviralEnrichment,
    #[serde(rename="culture-identification")]
    CultureIdentification
}
impl std::fmt::Display for Pipeline {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::PathogenDetection => write!(f, "pathogen-detection"),
            Self::PanviralEnrichment => write!(f, "panviral-enrichment"),
            Self::CultureIdentification => write!(f, "culture-identification"),

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
    pub pipeline: Pipeline,
    pub stage: String
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
            stage: schema.stage.clone()
        }
    }
}

