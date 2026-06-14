//! Shared worker context (S3-1).
//!
//! Built **once** at start-up and shared (`Arc`) across every job runner. Holds the
//! resolved [`WorkerConfig`] and the Cerebro API + FS clients.
//!
//! Both clients are **blocking** (`reqwest::blocking`); a job runner must therefore
//! invoke them via [`tokio::task::spawn_blocking`] (cloning the client, which is
//! `Clone`, into the closure) so the async executor is never blocked — see
//! [`WorkerContext::api`] / [`WorkerContext::fs`].

use std::sync::Arc;

use cerebro_client::client::CerebroClient;
use cerebro_fs::client::FileSystemClient;

use crate::config::WorkerConfig;
use crate::error::WorkerError;

pub struct WorkerContext {
    pub config: WorkerConfig,
    api: Option<CerebroClient>,
    fs: Option<FileSystemClient>,
}

impl WorkerContext {
    /// Build the context from config, tolerating an unconfigured API/FS so the
    /// worker still boots (for health/ping smoke tests). Lifecycle runners surface
    /// [`WorkerError::ClientNotConfigured`] until the `CEREBRO_*` env is supplied.
    ///
    /// Synchronous on purpose: call it from `spawn_blocking` so the blocking
    /// `reqwest` clients are constructed off the async executor.
    pub fn build(config: WorkerConfig) -> Arc<Self> {
        let (api, fs) = match &config.api_url {
            Some(url) => match CerebroClient::new(
                url,
                config.api_token.clone(),
                false,
                config.danger_invalid_certificate,
                config.api_token_file.clone(),
                config.team.clone(),
                config.db.clone(),
                config.project.clone(),
            ) {
                Ok(api) => {
                    let fs = FileSystemClient::with_config(&api, config.fs_config());
                    tracing::info!(api_url = %url, "Cerebro API + FS clients ready");
                    (Some(api), Some(fs))
                }
                Err(e) => {
                    tracing::warn!(error = %e, "failed to build Cerebro API client; lifecycle runners will error until configured");
                    (None, None)
                }
            },
            None => {
                tracing::warn!("CEREBRO_API_URL not set; lifecycle runners will error until configured");
                (None, None)
            }
        };

        Arc::new(Self { config, api, fs })
    }

    /// The API client (lifecycle endpoints). Clone the returned client into a
    /// `spawn_blocking` closure to call it — it is blocking.
    pub fn api(&self) -> Result<&CerebroClient, WorkerError> {
        self.api.as_ref().ok_or(WorkerError::ClientNotConfigured("API"))
    }

    /// The FS client (capture / verify / restore / physical storage ops). Clone the
    /// returned client into a `spawn_blocking` closure to call it — it is blocking.
    pub fn fs(&self) -> Result<&FileSystemClient, WorkerError> {
        self.fs.as_ref().ok_or(WorkerError::ClientNotConfigured("FS"))
    }
}
