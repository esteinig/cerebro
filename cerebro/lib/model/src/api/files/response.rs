use serde::Deserialize;
use crate::api::files::model::SeaweedFile;

#[derive(Deserialize)]
pub struct RegisterFileResponse {
    pub status: String,
    pub message: String,
    pub data: RegisterFileResponseData
}
#[derive(Deserialize)]
pub struct RegisterFileResponseData {
}

#[derive(Deserialize)]
pub struct ListFilesResponse {
    pub status: String,
    pub message: String,
    pub data: Vec<SeaweedFile>
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
