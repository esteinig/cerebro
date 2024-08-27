use std::path::PathBuf;

use reqwest::StatusCode;
use thiserror::Error;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum TowerError {
    #[error(transparent)]
    CerebroModelError(#[from] cerebro_model::api::cerebro::model::ModelError),

    #[error(transparent)]
    IoError(#[from] std::io::Error),
    
    #[error(transparent)]
    SerdeFailure(#[from] serde_json::Error),

    #[error("Failed to execute system process: {0}")]
    TowerProcessFailed(String),

    #[error("Reqwest error: {0}")]
    ReqwestFailure(#[from] reqwest::Error),

    #[error("Token file not found")]
    TokenFileNotFound,

    #[error("Token is missing")]
    TokenMissing,

    #[error("Failed to configure team")]
    RequireTeamNotConfigured,

    #[error("Data response failure: {1}")]
    DataResponseFailure(StatusCode, String),

    #[error("Failed response: {0}")]
    ResponseFailure(StatusCode),
}
