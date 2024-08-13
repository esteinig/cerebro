use serde::{Deserialize, Serialize};

use super::model::StagedSample;

#[derive(Serialize, Deserialize)]
pub struct RegisterStagedSampleResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<String>>
}
impl RegisterStagedSampleResponse {
    pub fn success(ids: Vec<String>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Sample staged successfully"),
            data: Some(ids)
        }
    }
    pub fn conflict(ids: Vec<String>) -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Staged sample with unique identifier already exists in database"),
            data: Some(ids)
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database operation: {}", error_message),
            data: None
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No staged samples entries found in database"),
            data: None
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct ListStagedSamplesResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<StagedSample>>
}
impl ListStagedSamplesResponse {
    pub fn success(samples: Vec<StagedSample>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Staged sample entries found in database"),
            data: Some(samples)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No staged samples entries found in database"),
            data: None
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct DeleteStagedSampleResponse {
    pub status: String,
    pub message: String,
    pub data: Option<StagedSample>
}
impl DeleteStagedSampleResponse {
    pub fn success(sample: StagedSample) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Staged sample entries deleted from database"),
            data: Some(sample)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No staged sample entries found in database"),
            data: None
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None
        }
    }
}