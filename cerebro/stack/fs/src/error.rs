use cerebro_client::error::HttpClientError;
use cerebro_model::api::cerebro::model::ModelError;
use reqwest::StatusCode;
use serde_json::Error as SerdeError;
use std::{io::Error as IoError, path::PathBuf};
use thiserror::Error;

/// Errors that can occur during oepration of SeaweedFS
#[derive(Error, Debug)]
pub enum WeedError {
    #[error(transparent)]
    ParseError(#[from] SerdeError),
    #[error(transparent)]
    IoError(#[from] IoError),
    #[error("received more than one upload response")]
    SingleUploadResponseError,
    /// Represents a failure to excecute a program in the weed pipeline
    #[error("Failed to execute '{0}' - is it installed?")]
    ProgramExecutionFailed(String),
    /// Represents a failure to excecute a process command in the weed pipeline
    #[error("Failed to run command, output is:\n{0}")]
    CommandExecutionFailed(String),
}

/// Errors that can occur during the download and installation of the `weed` executable
#[derive(Error, Debug)]
pub enum WeedDownloadError {
    #[error("failed to download the weed binary")]
    DownloadError,
    #[error("failed to extract the weed binary")]
    ExtractError,
    #[error("failed to set permissions for the weed binary")]
    PermissionError,
    #[error("failed to move the weed binary to the desired path")]
    MoveError,
}

/// Errors raised by the SeaweedFS filer HTTP client ([`crate::filer`]).
#[derive(Error, Debug)]
pub enum FilerError {
    /// Underlying HTTP/transport error.
    #[error("filer request failed")]
    Http(#[from] reqwest::Error),
    /// I/O error reading the local file or writing the downloaded object.
    #[error(transparent)]
    Io(#[from] IoError),
    /// Failed to deserialise the filer response body.
    #[error("failed to parse filer response")]
    Parse(#[from] SerdeError),
    /// The local file to upload does not exist.
    #[error("local file does not exist: {0}")]
    LocalFileMissing(String),
    /// The requested remote object was not found.
    #[error("remote object not found: {0}")]
    NotFound(String),
    /// The filer reported an application-level error during upload.
    #[error("filer upload error: {0}")]
    Upload(String),
    /// The filer returned an unexpected HTTP status.
    #[error("unexpected filer response status: {0}")]
    UnexpectedStatus(StatusCode),
}

#[derive(Error, Debug)]
pub enum FileSystemError {
    #[error("network error occurred")]
    NetworkError(#[from] reqwest::Error),
    #[error("cluster is unhealthy")]
    UnhealthyCluster,
    #[error("failed to extract file name")]
    FileNameExtraction,
    #[error("file does not exist: {0}")]
    FileDoesNotExist(PathBuf),
    #[error("failed to serialise run manifest: {0}")]
    ManifestSerialization(String),
    #[error("unexpected response status: {0}")]
    UnexpectedResponseStatus(StatusCode),
    #[error("I/O error occurred")]
    IoError(#[from] IoError),
    #[error("error during file download")]
    WeedError(#[from] WeedError),
    #[error("error during executable download")]
    WeedDownloadError(#[from] WeedDownloadError),
    #[error("error during file request")]
    HttpClientError(#[from] HttpClientError),
    #[error("error during model parsing")]
    ModelError(#[from] ModelError),
    #[error(transparent)]
    CsvError(#[from] csv::Error),
    /// Error during a SeaweedFS filer operation.
    #[error("error during filer operation")]
    FilerError(#[from] FilerError),
    /// A downloaded file failed BLAKE3 integrity verification.
    #[error("integrity verification failed for {path}: expected {expected}, got {actual}")]
    IntegrityMismatch {
        /// Local path of the file whose hash did not match.
        path: String,
        /// Hash registered with the Cerebro API.
        expected: String,
        /// Hash recomputed from the downloaded bytes.
        actual: String,
    },
    /// A download was requested without sufficient query parameters.
    #[error("invalid download query: provide --fid, or --run-id (optionally with --sample-id)")]
    InvalidDownloadQuery,
}
