use serde::{Deserialize, Serialize};

use super::model::ProductionPipeline;

#[derive(Serialize, Deserialize)]
pub struct RegisterPipelineResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>
}
impl RegisterPipelineResponse {
    pub fn success(id: &str) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Pipeline registered successfully"),
            data: Some(id.to_string())
        }
    }
    pub fn conflict(id: &str, name: &str, location: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("Pipeline with identifier '{id}' or name '{name}' and location '{location}' already exists in database"),
            data: Some(id.to_string())
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database operation: {}", error_message),
            data: None
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct ListPipelinesResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<ProductionPipeline>>
}
impl ListPipelinesResponse {
    pub fn success(pipelines: Vec<ProductionPipeline>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Pipeline registrations found in database"),
            data: Some(pipelines)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No pipeline registrations found in database"),
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
pub struct DeletePipelineResponse {
    pub status: String,
    pub message: String,
    pub data: Option<ProductionPipeline>
}
impl DeletePipelineResponse {
    pub fn success(pipeline: ProductionPipeline) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Pipeline registrations deleted from database"),
            data: Some(pipeline)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No pipeline registrations found in database"),
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
pub struct PingPipelineResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>
}
impl PingPipelineResponse {
    pub fn success(last_ping: String) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Pipeline last ping datetime updated"),
            data: Some(last_ping)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No pipeline registrations found in database"),
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