use serde::{Deserialize, Serialize};
use crate::api::files::model::SeaweedFile;

use super::model::FileTag;

#[derive(Serialize, Deserialize)]
pub struct RegisterFileResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>
}
impl RegisterFileResponse {
    pub fn success(id: &str) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File registered successfully"),
            data: Some(id.to_string())
        }
    }
    pub fn conflict(id: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("File with unique identifier already exists in database"),
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
pub struct ListFilesResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<SeaweedFile>>
}
impl ListFilesResponse {
    pub fn success(files: Vec<SeaweedFile>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File entries found in database"),
            data: Some(files)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No file entries found in database"),
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
pub struct UpdateTagsResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<FileTag>>
}
impl UpdateTagsResponse {
    pub fn success() -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Tags updated for files"),
            data: None
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Tags updated for files"),
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

/// Represents the output of a successful `weed` command execution.
#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WeedUploadResponse {
    pub file_name: String,
    pub url: String,
    pub fid: String,
    pub size: u64,
}


#[derive(Serialize, Deserialize)]
pub struct DeleteFileResponse {
    pub status: String,
    pub message: String,
    pub data: Option<SeaweedFile>
}
impl DeleteFileResponse {
    pub fn success(file: SeaweedFile) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File entries deleted from database"),
            data: Some(file)
        }
    }
    pub fn all_deleted() -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File entries deleted from database"),
            data: None
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No file entries found in database"),
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