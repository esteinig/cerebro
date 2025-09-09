use serde::{Deserialize, Serialize};
use crate::api::{cerebro::schema::PrefetchData, training::model::TrainingPrefetchRecord};

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
    pub identifier: String,
    pub name: String,
    pub prefetch: PrefetchData,
}

impl TrainingPrefetchData {
    pub fn from_query(record: TrainingPrefetchRecord, prefetch: PrefetchData) -> Self {
        Self {
            id: record.id,
            collection: record.collection,
            identifier: record.identifier,
            name: record.name,
            prefetch: prefetch
        }
    }
}