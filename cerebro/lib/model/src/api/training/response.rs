use serde::{Deserialize, Serialize};
use crate::api::{cerebro::schema::PrefetchData, training::{model::{TrainingPrefetchRecord, TrainingSessionRecord}, schema::TrainingRecord}};

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
    /// Timestamp when training started
    pub started: String,
    /// Timestamp when training completed (if any)
    pub completed: Option<String>,
    /// Number of records in the session
    pub records: Vec<TrainingRecord>,
}

impl TrainingSessionData {
    pub fn from_query(session: TrainingSessionRecord) -> Self {
        Self {
            id: session.id.clone(),
            user_name: session.user_name.clone(),
            user_id: session.user_id.clone(),
            started: session.started.clone(),
            completed: session.completed.clone(),
            records: session.records.clone(),
        }
    }
}