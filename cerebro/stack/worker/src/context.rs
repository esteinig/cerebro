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

/// The authenticated client pair, rebuilt atomically on (re-)login (S3-5 #3).
struct Clients {
    api: CerebroClient,
    fs: FileSystemClient,
}

pub struct WorkerContext {
    pub config: WorkerConfig,
    /// API + FS clients behind a lock so a periodic Bot re-login can rebuild them
    /// in place (decision #3) without a shared mutable token cell.
    clients: std::sync::RwLock<Option<Clients>>,
}

impl WorkerContext {
    /// Build the context from config, tolerating an unconfigured API/FS so the
    /// worker still boots (for health/ping smoke tests). Lifecycle runners surface
    /// [`WorkerError::ClientNotConfigured`] until the `CEREBRO_*` env is supplied.
    ///
    /// No login happens here — call [`WorkerContext::login`] after building (from an
    /// async context) to authenticate as the service Bot.
    pub fn build(config: WorkerConfig) -> Arc<Self> {
        let clients = Self::build_clients(&config);
        Arc::new(Self {
            config,
            clients: std::sync::RwLock::new(clients),
        })
    }

    /// Construct the API + FS clients from config (no login). `None` when
    /// `CEREBRO_API_URL` is unset or the client cannot be built.
    fn build_clients(config: &WorkerConfig) -> Option<Clients> {
        let url = config.api_url.as_ref()?;
        match CerebroClient::new(
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
                Some(Clients { api, fs })
            }
            Err(e) => {
                tracing::warn!(error = %e, "failed to build Cerebro API client; lifecycle runners will error until configured");
                None
            }
        }
    }

    /// The API client (lifecycle endpoints). Returns a clone (cheap) taken under a
    /// read lock so a concurrent re-login can swap the underlying client. Move the
    /// returned client into a `spawn_blocking` closure to call it — it is blocking.
    pub fn api(&self) -> Result<CerebroClient, WorkerError> {
        self.clients
            .read()
            .unwrap()
            .as_ref()
            .map(|c| c.api.clone())
            .ok_or(WorkerError::ClientNotConfigured("API"))
    }

    /// The FS client (capture / verify / restore / physical storage ops). Returns a
    /// clone taken under a read lock. Move it into a `spawn_blocking` closure — it
    /// is blocking.
    pub fn fs(&self) -> Result<FileSystemClient, WorkerError> {
        self.clients
            .read()
            .unwrap()
            .as_ref()
            .map(|c| c.fs.clone())
            .ok_or(WorkerError::ClientNotConfigured("FS"))
    }

    /// Perform the initial service-Bot login (S3-5 #5). Equivalent to
    /// [`relogin`](Self::relogin); a no-op when bot credentials aren't configured
    /// (the worker then runs in static-token or unauthenticated mode).
    pub async fn login(self: &Arc<Self>) -> Result<(), WorkerError> {
        self.relogin().await
    }

    /// Re-authenticate as the service Bot and rebuild the API + FS clients in place
    /// (decision #3). The fresh token is written to `api_token_file` by the client
    /// and the rebuilt clients are swapped under the write lock, so all subsequent
    /// `api()`/`fs()` clones use the new token. Serialised by the lock so concurrent
    /// callers don't stampede the login endpoint. No-op without bot credentials.
    pub async fn relogin(self: &Arc<Self>) -> Result<(), WorkerError> {
        let (Some(email), Some(password)) = (
            self.config.bot_email.clone(),
            self.config.bot_password.clone(),
        ) else {
            return Ok(());
        };
        let config = self.config.clone();
        let clients = self
            .run_blocking(move || {
                let url = config
                    .api_url
                    .as_ref()
                    .ok_or(WorkerError::ClientNotConfigured("API"))?;
                let api = CerebroClient::new(
                    url,
                    config.api_token.clone(),
                    false,
                    config.danger_invalid_certificate,
                    config.api_token_file.clone(),
                    config.team.clone(),
                    config.db.clone(),
                    config.project.clone(),
                )
                .map_err(|e| WorkerError::Api(e.to_string()))?;
                // Bot login: ?role=Bot yields the long-lived, role-scoped token.
                api.login_user(&email, Some(password), true)
                    .map_err(|e| WorkerError::Api(format!("bot login failed: {e}")))?;
                let fs = FileSystemClient::with_config(&api, config.fs_config());
                Ok(Clients { api, fs })
            })
            .await?;
        *self.clients.write().unwrap() = Some(clients);
        tracing::info!("authenticated as service bot; API + FS clients (re)built");
        Ok(())
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
        self.enqueue_at(kind, args, queue, None, reserve_for, None)
            .await
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
