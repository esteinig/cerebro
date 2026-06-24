//! CIQA API response envelope + read DTOs. Mirrors `training::response::TrainingResponse`.

use serde::{Deserialize, Serialize};

use crate::api::ciqa::model::{QualityControlBaseline, QualityControlDatasetRecord};

/// Uniform CIQA API envelope (mirror of `TrainingResponse<T>`).
#[derive(Debug, Serialize, Deserialize)]
pub struct CiqaResponse<T> {
    pub status: String,
    pub message: String,
    pub data: Option<T>,
}

impl<T> CiqaResponse<T> {
    pub fn ok(data: T) -> Self {
        Self { status: "ok".into(), message: "ok".into(), data: Some(data) }
    }
    pub fn created(message: &str) -> Self {
        Self { status: "ok".into(), message: message.into(), data: None }
    }
    pub fn error(message: &str) -> Self {
        Self { status: "fail".into(), message: message.into(), data: None }
    }
    pub fn not_found(message: &str) -> Self {
        Self { status: "fail".into(), message: message.into(), data: None }
    }
    pub fn conflict(message: &str) -> Self {
        Self { status: "fail".into(), message: message.into(), data: None }
    }
}

/// Listing DTO for datasets (metadata only; the `PrefetchData` stays in GridFS).
#[derive(Debug, Serialize, Deserialize)]
pub struct CiqaDatasetList {
    pub datasets: Vec<QualityControlDatasetRecord>,
}

/// Listing DTO for baselines.
#[derive(Debug, Serialize, Deserialize)]
pub struct CiqaBaselineList {
    pub baselines: Vec<QualityControlBaseline>,
}
