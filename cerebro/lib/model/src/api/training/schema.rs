use std::path::Path;

use serde::{Deserialize, Serialize};
use uuid::Uuid;
use crate::api::cerebro::{model::ModelError, schema::{PrefetchData, SampleType, TestResult}};

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
pub struct TrainingPrefetchOverview {
    pub collection: String,
    pub description: String,
    pub samples: i64,
    pub session_id: Option<String>, // last completed session for this user in this collection
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingRecord {
    pub id: String, 
    pub data_id: String,    // matches identifier of TrainingPrefetchData
    pub result: TestResult,
    pub sample_name: Option<String>,
    pub sample_type: Option<SampleType>,
    pub candidates: Option<Vec<String>>,
    pub reference_result: Option<TestResult>,
    pub reference_candidates: Option<Vec<String>>,
    pub exclude_lod: Option<bool>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateTrainingSession {
    pub collection: String,
    pub shuffle: bool
}

#[derive(Debug, Clone, Deserialize)]
pub struct PatchTrainingRecord {
    /// TrainingSession.id to update
    pub session_id: String,
    /// TrainingRecord.id to update
    pub record_id: String,
    /// New test result
    pub result: TestResult,
    /// Optional candidates list
    pub candidates: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryTrainingPrefetch {
    /// Optional filter by collection
    pub collection: Option<String>,
    /// Optional filter by identifier
    pub id: Option<String>,
    /// Optional filter by name
    pub name: Option<String>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryTrainingData {
    /// Record index in the current session
    pub record: u64,
}