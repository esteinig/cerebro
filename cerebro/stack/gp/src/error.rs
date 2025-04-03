use std::path::PathBuf;

use thiserror::Error;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum GptError {
    #[error(transparent)]
    SerdeJsonError(#[from] serde_json::Error),
    #[error(transparent)]
    IoError(#[from] std::io::Error),
    #[error(transparent)]
    CerebroModelError(#[from] cerebro_model::api::cerebro::model::ModelError),
    #[error(transparent)]
    NifflerError(#[from] niffler::Error),
    #[error(transparent)]
    CsvError(#[from] csv::Error)
}
