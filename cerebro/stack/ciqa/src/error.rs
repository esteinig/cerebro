use std::path::PathBuf;

use thiserror::Error;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum CiqaError {
    #[error(transparent)]
    IoError(#[from] std::io::Error),
    #[error(transparent)]
    SerdeError(#[from] serde_json::Error),
    #[error(transparent)]
    CerebroClientError(#[from] cerebro_client::error::HttpClientError),
    #[error(transparent)]
    NvlmError(#[from] nvml_wrapper::error::NvmlError),
    #[error(transparent)]
    CerebroModelError(#[from] cerebro_model::api::cerebro::model::ModelError),
    #[error("plotters crate error: {0}")]
    PlottersError(#[from] Box<dyn std::error::Error + Send + Sync>), 
    #[error(transparent)]
    NifflerError(#[from] niffler::Error),
    #[error(transparent)]
    CsvError(#[from] csv::Error),
    #[error(transparent)]
    MetaGptError(#[from] meta_gpt::error::GptError),
    #[error("failed to extract sample identifier from filename (format: sample_id.model.json)")]
    SampleIdentifierNotFound,
    #[error("failed to extract model name from filename (format: sample_id.model.json)")]
    ModelNameNotFound,
    #[error("failed to extract file stem ({0})")]
    FileStemNotFound(PathBuf),
    #[error("failed to convert OsString to String: {0}")]
    FileNameConversionError(PathBuf),
    #[error("error during strip-plot creation: {0}")]
    StripPlotError(String),
    #[error("no columns found for plate plot")]
    NoColumnsFound,
    #[error("invalid hex color string: {0}")]
    InvalidHexColor(String),
    #[error("invalid diagnostic outcome token in review matrix: '{0}' (expected one of TP, TN, FP, FN, NA, Indeterminate, Control, Unknown)")]
    InvalidOutcomeToken(String),
    #[error("review matrix file is empty or missing a header row: {0}")]
    EmptyReviewMatrix(PathBuf),
    #[error("review matrix must contain a sample column plus at least one reviewer column: {0}")]
    NoReviewerColumns(PathBuf),
    #[error("sample identifier must be specified when not using prefetch data")]
    SampleIdentifierMissing, 
}

impl<T> From<plotters::drawing::DrawingAreaErrorKind<T>> for CiqaError
where
    T: std::error::Error + Send + Sync + 'static,
{
    fn from(err: plotters::drawing::DrawingAreaErrorKind<T>) -> Self {
        CiqaError::PlottersError(Box::new(err))
    }
}