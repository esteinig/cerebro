//! Worker configuration, read from the environment (S3-1).
//!
//! The worker is a standalone process; everything it needs to reach the Cerebro
//! API, Cerebro FS and Faktory comes from `CEREBRO_*` / `FAKTORY_URL` env vars so
//! it can be configured purely through the deployment (docker secrets / env_file),
//! mirroring how the server and tower are configured.

use cerebro_fs::config::{FsAccessMode, FsConfig};

/// Resolved worker configuration.
#[derive(Debug, Clone)]
pub struct WorkerConfig {
    // --- Cerebro API (lifecycle endpoints) ---
    pub api_url: Option<String>,
    pub api_token: Option<String>,
    pub api_token_file: Option<std::path::PathBuf>,
    pub team: Option<String>,
    pub db: Option<String>,
    pub project: Option<String>,
    pub danger_invalid_certificate: bool,

    // --- Cerebro FS (physical storage ops) ---
    pub fs_master_url: String,
    pub fs_master_port: String,
    pub fs_filer_url: String,
    pub fs_access: FsAccessMode,

    // --- Worker runtime ---
    /// Faktory queues this worker consumes, in priority order.
    pub queues: Vec<String>,
    /// `host:port` for the worker's health/metrics HTTP server.
    pub metrics_addr: String,
    /// Run the deep integrity gate (download + BLAKE3) before committing a tier
    /// move (S3-2a). Heavy for large artefacts; off by default — deep verification
    /// is the scheduled verify worker's job (S3-3a). Per-job override via the
    /// `tier_move` arg `verify: true`.
    pub verify_on_move: bool,
    /// Dev/test simulation for the restore executor (S3-3b): when set, an archival
    /// restore is treated as ready this many seconds after it was requested,
    /// letting the state machine be exercised end-to-end without a real S3 Glacier
    /// integration. Unset in production (where the real provider drives readiness).
    pub restore_simulate_seconds: Option<i64>,
}

fn env_opt(key: &str) -> Option<String> {
    std::env::var(key).ok().filter(|v| !v.is_empty())
}

fn env_or(key: &str, default: &str) -> String {
    env_opt(key).unwrap_or_else(|| default.to_string())
}

fn env_bool(key: &str) -> bool {
    matches!(env_opt(key).as_deref(), Some("true") | Some("1") | Some("yes"))
}

impl WorkerConfig {
    /// Read configuration from the environment, applying sane defaults.
    pub fn from_env() -> Self {
        let access = match env_opt("CEREBRO_FS_ACCESS").as_deref() {
            Some("filer") => FsAccessMode::Filer,
            Some("weed") => FsAccessMode::Weed,
            _ => FsAccessMode::default(),
        };

        let queues = env_or("CEREBRO_WORKER_QUEUES", "default,lifecycle,maintenance")
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect::<Vec<_>>();

        Self {
            api_url: env_opt("CEREBRO_API_URL"),
            api_token: env_opt("CEREBRO_API_TOKEN"),
            api_token_file: env_opt("CEREBRO_API_TOKEN_FILE").map(Into::into),
            team: env_opt("CEREBRO_TEAM"),
            db: env_opt("CEREBRO_DB"),
            project: env_opt("CEREBRO_PROJECT"),
            danger_invalid_certificate: env_bool("CEREBRO_DANGER_INVALID_CERTIFICATE"),

            fs_master_url: env_or("CEREBRO_FS_URL", "http://localhost"),
            fs_master_port: env_or("CEREBRO_FS_PORT", "9333"),
            fs_filer_url: env_or("CEREBRO_FS_FILER_URL", "http://localhost:8888"),
            fs_access: access,

            queues: if queues.is_empty() { vec!["default".to_string()] } else { queues },
            metrics_addr: env_or("CEREBRO_WORKER_METRICS_ADDR", "0.0.0.0:9464"),
            verify_on_move: env_bool("CEREBRO_WORKER_VERIFY_ON_MOVE"),
            restore_simulate_seconds: env_opt("CEREBRO_RESTORE_SIMULATE_SECONDS")
                .and_then(|v| v.parse::<i64>().ok()),
        }
    }

    /// Build the [`FsConfig`] for the FS client from this config.
    pub fn fs_config(&self) -> FsConfig {
        FsConfig {
            master_url: self.fs_master_url.clone(),
            master_port: self.fs_master_port.clone(),
            localhost: true,
            filer_url: self.fs_filer_url.clone(),
            access: self.fs_access.clone(),
            danger_invalid_certificate: self.danger_invalid_certificate,
        }
    }
}
