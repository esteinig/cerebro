use std::{fs::File, io::BufReader, path::Path};

use serde::{Deserialize, Serialize};
use crate::api::{cerebro::{model::ModelError, schema::PrefetchData}, training::{model::{TrainingPrefetchRecord, TrainingResult, TrainingSessionRecord}, schema::TrainingRecord}};

#[derive(Debug, Serialize, Deserialize)]
pub struct TrainingResponse<T> {
    pub status: String,
    pub message: String,
    pub data: Option<T>,
}

impl<T> TrainingResponse<T> {
    pub fn ok(data: T) -> Self {
        Self { status: "ok".to_string(), message: "ok".into(), data: Some(data) }
    }
    pub fn completed() -> Self {
        Self { status: "ok".to_string(), message: "ok".into(), data: None }
    }
    pub fn created(message: &str) -> Self {
        Self { status: "ok".to_string(), message: message.into(), data: None }
    }
    pub fn error(message: &str) -> Self {
        Self { status: "fail".to_string(), message: message.into(), data: None }
    }
    pub fn not_found(message: &str) -> Self {
        Self { status: "fail".to_string(), message: message.into(), data: None }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TrainingPrefetchData {
    pub id: String,
    pub collection: String,
    pub description: String,
    pub name: String,
    pub prefetch: PrefetchData,
}

impl TrainingPrefetchData {
    pub fn from_query(record: TrainingPrefetchRecord, prefetch: PrefetchData) -> Self {
        Self {
            id: record.id,
            collection: record.collection,
            description: record.description,
            name: record.name,
            prefetch: prefetch
        }
    }
}


#[derive(Debug, Serialize, Deserialize)]
pub struct TrainingSessionData {
    /// Unique identifier for the training session
    pub id: String,
    /// User display name
    pub user_name: String,
    /// User unique identifier
    pub user_id: String,
    /// Training session 
    pub collection: String,
    /// Timestamp when training started
    pub started: String,
    /// Timestamp when training completed (if any)
    pub completed: Option<String>,
    /// Training result data
    pub result: Option<TrainingResult>,
    /// Number of records in the session
    pub records: Vec<TrainingRecord>,
}

impl TrainingSessionData {
    pub fn from_query(session: TrainingSessionRecord, evaluate: bool) -> Self {
        Self {
            id: session.id.clone(),
            user_name: session.user_name.clone(),
            user_id: session.user_id.clone(),
            started: session.started.clone(),
            collection: session.collection.clone(),
            completed: session.completed.clone(),
            result: if evaluate && session.completed.is_some() { Some(session.evaluate()) } else { None },
            records: session.records.clone(),
        }
    }
    pub fn to_json<P: AsRef<Path>>(&self, path: P) -> Result<(), ModelError> {
        let data = serde_json::to_string_pretty(&self).map_err(|err| ModelError::JsonSerialization(err))?;
        std::fs::write(path, data)?;
        Ok(())
    }
    pub fn from_json<P: AsRef<Path>>(path: P) -> Result<Self, ModelError> {
        let rdr = BufReader::new(File::open(&path)?);
        let config: Self = serde_json::from_reader(rdr).map_err(|err| ModelError::JsonDeserialization(err))?;
        Ok(config)
    }
}