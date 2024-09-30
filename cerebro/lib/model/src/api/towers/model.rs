use serde::{Serialize, Deserialize};
use super::schema::RegisterTowerSchema;

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
pub struct ProductionTower {
    pub id: String,
    pub date: String,
    pub name: String,
    pub location: String,
    pub last_ping: String,
    pub pipelines: Vec<Pipeline>,
    pub stage: String
}
impl ProductionTower {
    pub fn from_schema(schema: &RegisterTowerSchema) -> Self {
        Self {
            id: schema.id.clone(),
            date: schema.date.clone(),
            name: schema.name.clone(),
            location: schema.location.clone(),
            last_ping: schema.last_ping.clone(),
            pipelines: schema.pipelines.clone(),
            stage: schema.stage.clone()
        }
    }
}

