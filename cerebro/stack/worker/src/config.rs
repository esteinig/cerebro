//! Worker configuration, read from the environment.
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
    /// Service Bot login email. When set with `bot_password`, the worker
    /// logs in as this Bot at startup (and periodically re-logs in) to obtain a
    /// long-lived, role-scoped token, instead of relying on a static token.
    pub bot_email: Option<String>,
    pub bot_password: Option<String>,
    /// How often to re-authenticate as the Bot and rebuild the clients in place
    /// (seconds). Kept well under the bot token TTL (`access_max_age_bot`, 30d by
    /// default) so the token never lapses in practice. Default 7 days.
    pub relogin_interval_secs: u64,
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
    /// move. Heavy for large artefacts; off by default — deep verification
    /// is the scheduled verify worker's job. Per-job override via the
    /// `tier_move` arg `verify: true`.
    pub verify_on_move: bool,
    /// Dev/test simulation for the restore executor: when set, an archival
    /// restore is treated as ready this many seconds after it was requested,
    /// letting the state machine be exercised end-to-end without a real S3 Glacier
    /// integration. Unset in production (where the real provider drives readiness).
    pub restore_simulate_seconds: Option<i64>,
    /// Catalogue backup settings. `Some` only when both a backup MongoDB
    /// URI and a store path are configured; otherwise the `catalogue_backup`
    /// runner reports "not configured" and does nothing.
    pub backup: Option<BackupSettings>,
    /// Cold archival object-store settings. `Some` enables real archival of
    /// files moved to the Cold tier; `None` keeps the prior label-only behaviour.
    pub archive: Option<ArchiveSettings>,
}

/// Settings for the scheduled catalogue/audit backup.
#[derive(Debug, Clone)]
pub struct BackupSettings {
    /// Read-only MongoDB connection string `mongodump` runs against.
    pub mongo_uri: String,
    /// Database to dump (the Cerebro catalogue DB).
    pub mongo_db: String,
    /// Filesystem object-store root (a mounted backup volume / NFS). The S3
    /// backend will add a URL form behind the same store trait.
    pub store_root: std::path::PathBuf,
    /// Key prefix under the store for catalogue backups.
    pub prefix: String,
    /// Most-recent backups to always retain.
    pub keep_last: usize,
    /// Age floor (days) before an older-than-`keep_last` backup may be pruned.
    pub min_age_days: i64,
}

fn env_opt(key: &str) -> Option<String> {
    std::env::var(key).ok().filter(|v| !v.is_empty())
}

/// Read a secret from `<KEY>` directly, or from the file named by `<KEY>_FILE`
/// (the Docker-secret convention). The file's whole contents are used, trimmed.
fn env_or_file(key: &str) -> Option<String> {
    if let Some(v) = env_opt(key) {
        return Some(v);
    }
    let path = env_opt(&format!("{key}_FILE"))?;
    match std::fs::read_to_string(&path) {
        Ok(s) => {
            let s = s.trim().to_string();
            if s.is_empty() {
                None
            } else {
                Some(s)
            }
        }
        Err(e) => {
            tracing::warn!(%path, "failed to read secret file for {key}: {e}");
            None
        }
    }
}

fn env_or(key: &str, default: &str) -> String {
    env_opt(key).unwrap_or_else(|| default.to_string())
}

fn env_bool(key: &str) -> bool {
    matches!(
        env_opt(key).as_deref(),
        Some("true") | Some("1") | Some("yes")
    )
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
            bot_email: env_opt("CEREBRO_API_BOT_EMAIL"),
            bot_password: env_or_file("CEREBRO_API_BOT_PASSWORD"),
            relogin_interval_secs: env_opt("CEREBRO_WORKER_RELOGIN_SECS")
                .and_then(|v| v.parse::<u64>().ok())
                .unwrap_or(7 * 24 * 3600),
            team: env_opt("CEREBRO_TEAM"),
            db: env_opt("CEREBRO_DB"),
            project: env_opt("CEREBRO_PROJECT"),
            danger_invalid_certificate: env_bool("CEREBRO_DANGER_INVALID_CERTIFICATE"),

            fs_master_url: env_or("CEREBRO_FS_URL", "http://localhost"),
            fs_master_port: env_or("CEREBRO_FS_PORT", "9333"),
            fs_filer_url: env_or("CEREBRO_FS_FILER_URL", "http://localhost:8888"),
            fs_access: access,

            queues: if queues.is_empty() {
                vec!["default".to_string()]
            } else {
                queues
            },
            metrics_addr: env_or("CEREBRO_WORKER_METRICS_ADDR", "0.0.0.0:9464"),
            verify_on_move: env_bool("CEREBRO_WORKER_VERIFY_ON_MOVE"),
            restore_simulate_seconds: env_opt("CEREBRO_RESTORE_SIMULATE_SECONDS")
                .and_then(|v| v.parse::<i64>().ok()),
            backup: BackupSettings::from_env(),
            archive: ArchiveSettings::from_env(),
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

impl BackupSettings {
    /// Read backup settings from the environment. Returns `None` — backup
    /// disabled — unless both a backup MongoDB URI and a store path are provided,
    /// so the feature is strictly opt-in.
    pub fn from_env() -> Option<Self> {
        let mongo_uri = env_or_file("CEREBRO_BACKUP_MONGO_URI")?;
        let store_root = env_opt("CEREBRO_BACKUP_STORE_PATH")?;
        Some(Self {
            mongo_uri,
            mongo_db: env_or("CEREBRO_BACKUP_MONGO_DB", ""),
            store_root: store_root.into(),
            prefix: env_or("CEREBRO_BACKUP_PREFIX", "catalogue"),
            keep_last: env_opt("CEREBRO_BACKUP_KEEP_LAST")
                .and_then(|v| v.parse().ok())
                .unwrap_or(14),
            min_age_days: env_opt("CEREBRO_BACKUP_MIN_AGE_DAYS")
                .and_then(|v| v.parse().ok())
                .unwrap_or(7),
        })
    }
}

/// Cold archival object-store settings.
#[derive(Debug, Clone)]
pub struct ArchiveSettings {
    /// Key prefix under the cold store for archived objects.
    pub prefix: String,
    /// Selected cold-store backend.
    pub backend: ArchiveBackend,
    /// Grace period (days) after archival before the redundant local copy is
    /// reclaimed. Anchored on `tier_moved_at`. From
    /// `CEREBRO_ARCHIVE_LOCAL_GRACE_DAYS` (default 7).
    pub local_grace_days: i64,
}

/// Cold-store backend: a local/NFS directory (the compile-safe default) or an
/// S3-compatible store (requires building with `--features s3`).
#[derive(Debug, Clone)]
pub enum ArchiveBackend {
    Filesystem {
        root: std::path::PathBuf,
    },
    S3 {
        endpoint: String,
        region: String,
        bucket: String,
        access_key: String,
        secret_key: String,
    },
}

impl ArchiveSettings {
    /// Read archival settings from the environment. Returns `None` —
    /// label-only tiering, the prior behaviour — unless a backend is configured.
    /// An S3 endpoint selects the S3 backend; otherwise a store path selects the
    /// filesystem backend (a mounted cold disk / NFS).
    pub fn from_env() -> Option<Self> {
        let prefix = env_or("CEREBRO_ARCHIVE_PREFIX", "archive");
        let local_grace_days = env_opt("CEREBRO_ARCHIVE_LOCAL_GRACE_DAYS")
            .and_then(|v| v.parse().ok())
            .unwrap_or(7);
        if let Some(endpoint) = env_opt("CEREBRO_ARCHIVE_S3_ENDPOINT") {
            let bucket = env_opt("CEREBRO_ARCHIVE_S3_BUCKET")?;
            return Some(Self {
                prefix,
                local_grace_days,
                backend: ArchiveBackend::S3 {
                    endpoint,
                    region: env_or("CEREBRO_ARCHIVE_S3_REGION", "us-east-1"),
                    bucket,
                    access_key: env_or_file("CEREBRO_ARCHIVE_S3_ACCESS_KEY").unwrap_or_default(),
                    secret_key: env_or_file("CEREBRO_ARCHIVE_S3_SECRET_KEY").unwrap_or_default(),
                },
            });
        }
        let root = env_opt("CEREBRO_ARCHIVE_STORE_PATH")?;
        Some(Self {
            prefix,
            local_grace_days,
            backend: ArchiveBackend::Filesystem { root: root.into() },
        })
    }

    /// Open the configured cold object store. The archive *key* already carries the
    /// prefix (see `archive::archive_key_for`), so the store itself is given no
    /// extra prefix.
    pub fn open_store(&self) -> anyhow::Result<Box<dyn crate::backup::ObjectStore>> {
        match &self.backend {
            ArchiveBackend::Filesystem { root } => Ok(Box::new(
                crate::backup::FilesystemObjectStore::new(root.clone()),
            )),
            ArchiveBackend::S3 {
                endpoint,
                region,
                bucket,
                access_key,
                secret_key,
            } => {
                #[cfg(feature = "s3")]
                {
                    let store = crate::backup::S3ObjectStore::new(
                        endpoint, region, bucket, access_key, secret_key, "",
                    )?;
                    Ok(Box::new(store))
                }
                #[cfg(not(feature = "s3"))]
                {
                    let _ = (endpoint, region, bucket, access_key, secret_key);
                    anyhow::bail!(
                        "S3 archive backend requires building cerebro-worker with --features s3"
                    )
                }
            }
        }
    }
}
