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

    /// Run a blocking client/storage call on the blocking thread pool, off the
    /// async executor. This is the S3-1 concurrency decision made concrete: the
    /// Faktory protocol stays async, but lifecycle work is blocking and offloaded
    /// here, so the executor is never parked on I/O/CPU/subprocess work.
    ///
    /// Clone the client out of the context first, then move it into `f`:
    ///
    /// ```ignore
    /// let api = ctx.api()?.clone();
    /// let file = ctx.run_blocking(move || {
    ///     api.get_file(&id).map_err(|e| WorkerError::Api(e.to_string()))
    /// }).await?;
    /// ```
    pub async fn run_blocking<T, F>(&self, f: F) -> Result<T, WorkerError>
    where
        F: FnOnce() -> Result<T, WorkerError> + Send + 'static,
        T: Send + 'static,
    {
        tokio::task::spawn_blocking(f)
            .await
            .map_err(|e| WorkerError::Other(format!("blocking task join error: {e}")))?
    }

    /// Enqueue a Faktory job (producer side). Used by `*_scan` runners that fan out
    /// one job per due item. Opens a short-lived Faktory client from `FAKTORY_URL`
    /// (the same connection the server's enqueue path uses). `reserve_for` bounds
    /// how long the job may run before Faktory re-delivers it — set generously for
    /// long-running moves.
    pub async fn enqueue(
        &self,
        kind: &str,
        args: serde_json::Value,
        queue: &str,
        reserve_for: Option<std::time::Duration>,
    ) -> Result<(), WorkerError> {
        self.enqueue_at(kind, args, queue, None, reserve_for, None).await
    }

    /// Enqueue a Faktory job, optionally scheduled for a future time (`at`) and with
    /// an explicit `retry` count. This is the mechanism behind **poll-by-re-enqueue**
    /// (S3-3b restore): instead of parking a worker on a multi-hour archival thaw,
    /// the job checks status and, if not ready, re-enqueues *itself* at
    /// `now + poll_interval`. Poll jobs set `retry = 0` so Faktory's own retry never
    /// races the explicit poll chain.
    pub async fn enqueue_at(
        &self,
        kind: &str,
        args: serde_json::Value,
        queue: &str,
        at: Option<chrono::DateTime<chrono::Utc>>,
        reserve_for: Option<std::time::Duration>,
        retry: Option<isize>,
    ) -> Result<(), WorkerError> {
        let mut client = faktory::Client::connect()
            .await
            .map_err(|e| WorkerError::Other(format!("faktory connect: {e}")))?;
        let mut job = faktory::Job::new(kind, vec![args]).on_queue(queue);
        job.at = at;
        job.reserve_for = reserve_for;
        if let Some(r) = retry {
            job.retry = Some(r);
        }
        client
            .enqueue(job)
            .await
            .map_err(|e| WorkerError::Other(format!("faktory enqueue: {e}")))?;
        Ok(())
    }
}
