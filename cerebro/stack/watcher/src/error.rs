use std::path::PathBuf;

use thiserror::Error;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum WatcherError {
    
    #[error(transparent)]
    CerebroClientError(#[from] cerebro_client::error::HttpClientError),
    #[error(transparent)]
    CerebroModelError(#[from] cerebro_model::api::cerebro::model::ModelError),
    #[error(transparent)]
    CerebroSlackError(#[from] cerebro_model::slack::SlackError),
    #[error(transparent)]
    FileSystemClientError(#[from] cerebro_fs::error::FileSystemError),
    
    #[error(transparent)]
    NotifyError(#[from] notify::Error),
    #[error(transparent)]
    IoError(#[from] std::io::Error),

    #[error("failed to extract pattern match for sample identifier from: {0}")]
    GlobMatchSampleIdentifier(PathBuf),
    #[error("failed to extract an entry from the globbed walk through directory: {0}")]
    GlobWalk(PathBuf),
    #[error("failed to create a glob matcher for pattern: {0}")]
    GlobCreate(String),
    #[error("failed to find paired files for sample: {0}")]
    GlobPairedFiles(String),

    #[error("failed to detect read directory: {0}")]
    InvalidDirectoryStructure(PathBuf),
    #[error("failed to detect latest 'Analysis' directory in: {0}")]
    InvalidLatestAnalysisDirectory(PathBuf),
    #[error("failed to detect latest 'Alignment' directory in: {0}")]
    InvalidLatestAlignmentDirectory(PathBuf),


    #[error("failed to create a watcher configuration - were all required arguments provided? (--name & --location & --format)")]
    InvalidWatcherConfigArgs,
    #[error("failed to create a watcher configuration - were all required arguments provided? (--id | --json)")]
    WatcherIdentifierArgNotFound
}
