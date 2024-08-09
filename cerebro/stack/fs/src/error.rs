use cerebro_client::error::HttpClientError;
use thiserror::Error;
use std::{io::Error as IoError, path::PathBuf};
use serde_json::Error as SerdeError;
use reqwest::StatusCode;

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

/// Errors that may occur while executing the `weed` upload command
#[derive(Error, Debug)]
pub enum WeedUploadError {
    #[error("received more than one upload response")]
    SingleUploadResponseError,
    #[error("failed to parse command output")]
    ParseError(#[from] SerdeError),
    #[error("I/O error occurred")]
    IoError(#[from] IoError),
    /// Represents a failure to excecute a program in the simulation pipeline
    #[error("Failed to execute '{0}' - is it installed?")]
    ProgramExecutionFailed(String),
    /// Represents a failure to excecute a process command in the simulation pipeline
    #[error("Failed to run command, output is:\n{0}")]
    CommandExecutionFailed(String),
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
    #[error("unexpected response status: {0}")]
    UnexpectedResponseStatus(StatusCode),
    #[error("I/O error occurred")]
    IoError(#[from] IoError),
    #[error("error during file download")]
    WeedDownloadError(#[from] WeedDownloadError),
    #[error("error during file upload")]
    WeedUploadError(#[from] WeedUploadError),
    #[error("error during file request")]
    HttpClientError(#[from] HttpClientError),
}
