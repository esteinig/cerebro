use std::path::Path;

use serde::{Deserialize, Serialize};
use uuid::Uuid;
use crate::api::cerebro::{model::ModelError, schema::{PrefetchData, TestResult}};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateTrainingPrefetch {
    pub id: String,
    pub name: String,
    pub collection: String,
    pub description: String,
    pub prefetch: PrefetchData,
}
impl CreateTrainingPrefetch {
    pub fn from_file(path: &Path, collection: &str, description: &str) -> Result<Self, ModelError> {

        let prefetch_data = PrefetchData::from_json(path)?;

        Ok(Self {
            id: Uuid::new_v4().to_string(),
            collection: collection.to_string(),
            description: description.to_string(),
            name: prefetch_data.config.sample.clone(),
            prefetch: prefetch_data
        })

    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TrainingResult {
    Positive,
    Negative
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingRecord {
    pub id: String,
    pub collection: String,
    pub result: TestResult,
    pub candidates: Option<Vec<String>>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateTrainingSession {
    pub id: String,
    pub user_name: String,
    pub user_id: String,
    pub started: String,
    pub completed: Option<String>,
    pub records: Vec<TrainingRecord>
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryTrainingPrefetch {
    /// Optional filter by collection
    pub collection: Option<String>,
    /// Optional filter by identifier
    pub identifier: Option<String>,
    /// Optional filter by name
    pub name: Option<String>
}
