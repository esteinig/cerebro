use plotters::prelude::DrawingBackend;
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
    CerebroModelError(#[from] cerebro_model::api::cerebro::model::ModelError),
    #[error("Plotters error: {0}")]
    PlottersError(#[from] Box<dyn std::error::Error + Send + Sync>), 
    #[error(transparent)]
    NifflerError(#[from] niffler::Error),
    #[error(transparent)]
    CsvError(#[from] csv::Error),
    #[error(transparent)]
    GptError(#[from] cerebro_gp::error::GptError),
}

impl<T> From<plotters::drawing::DrawingAreaErrorKind<T>> for CiqaError
where
    T: std::error::Error + Send + Sync + 'static,
{
    fn from(err: plotters::drawing::DrawingAreaErrorKind<T>) -> Self {
        CiqaError::PlottersError(Box::new(err))
    }
}