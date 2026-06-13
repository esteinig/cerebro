//! Storage-layer configuration for `cerebro-fs`.
//!
//! [`FsConfig`] centralises the endpoints and access policy used to reach a
//! SeaweedFS deployment. It replaces the previous ad-hoc `(fs_url, fs_port,
//! localhost)` triplet that was threaded through the client, giving every
//! consumer (`cerebro-fs`, `cerebro-watcher`, and `cerebro-tower`)
//! a single, environment-driven source of truth for the master, filer, and
//! (later) S3 endpoints.
//!
//! # Access modes
//!
//! Object I/O can be performed two ways, selected by [`FsAccessMode`]:
//!
//! * [`FsAccessMode::Weed`] — shell out to the `weed` binary (master/volume,
//!   fid-addressed). This is the default. The `weed` client transparently
//!   chunks very large files (via its `-maxMB` manifest), so multi-gigabyte
//!   sequencing read sets upload safely without any in-process buffering.
//! * [`FsAccessMode::Filer`] — talk to the SeaweedFS Filer HTTP API
//!   (path-addressed; see [`crate::filer`]). Uploads/downloads stream from/to
//!   disk and the filer auto-chunks large objects, so this path is also safe
//!   for large files. In the Cerebro deployment the filer is backed by MongoDB
//!   (reusing the existing `cerebro-database`); see the deployment notes.

use std::env;

/// Default master base URL when none is supplied via flag or environment.
pub const DEFAULT_MASTER_URL: &str = "http://localhost";
/// Default master port (SeaweedFS master HTTP/health port).
pub const DEFAULT_MASTER_PORT: &str = "9333";
/// Default filer base URL (SeaweedFS filer HTTP API).
pub const DEFAULT_FILER_URL: &str = "http://localhost:8888";

/// Object I/O access mode for the storage client.
///
/// The variant names map onto the lower-cased values accepted on the command
/// line and in the `CEREBRO_FS_ACCESS` environment variable (`weed`, `filer`).
#[derive(Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum FsAccessMode {
    /// Shell out to the `weed` binary (fid-addressed, master/volume).
    Weed,
    /// SeaweedFS Filer HTTP API (path-addressed, streaming, MongoDB-backed store).
    Filer,
}

impl Default for FsAccessMode {
    fn default() -> Self {
        FsAccessMode::Weed
    }
}

impl std::fmt::Display for FsAccessMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FsAccessMode::Weed => write!(f, "weed"),
            FsAccessMode::Filer => write!(f, "filer"),
        }
    }
}

/// Connection configuration for a SeaweedFS deployment.
#[derive(Clone, Debug)]
pub struct FsConfig {
    /// Base URL of the SeaweedFS master (scheme + host), e.g. `http://localhost`.
    pub master_url: String,
    /// Master port, appended to `master_url` when [`FsConfig::localhost`] is set.
    pub master_port: String,
    /// Whether the master is addressed as `host:port` (localhost/dev) rather than
    /// a routed URL behind a reverse proxy (production).
    pub localhost: bool,
    /// Base URL of the SeaweedFS Filer HTTP API, e.g. `http://localhost:8888`.
    pub filer_url: String,
    /// Object I/O access mode.
    pub access: FsAccessMode,
    /// Accept invalid TLS certificates. Intended for development only.
    pub danger_invalid_certificate: bool,
}

impl Default for FsConfig {
    fn default() -> Self {
        Self {
            master_url: DEFAULT_MASTER_URL.to_string(),
            master_port: DEFAULT_MASTER_PORT.to_string(),
            localhost: true,
            filer_url: DEFAULT_FILER_URL.to_string(),
            access: FsAccessMode::default(),
            danger_invalid_certificate: false,
        }
    }
}

impl FsConfig {
    /// Construct a configuration for the `weed` (master/volume) access mode.
    ///
    /// This mirrors the historical behaviour of `FileSystemClient::new` and is
    /// used by the backwards-compatible constructor so existing callers keep
    /// working unchanged.
    ///
    /// # Examples
    ///
    /// ```
    /// use cerebro_fs::config::{FsConfig, FsAccessMode};
    ///
    /// let cfg = FsConfig::weed("http://localhost", "9333", true);
    /// assert_eq!(cfg.access, FsAccessMode::Weed);
    /// assert_eq!(cfg.master_health_url(), "http://localhost:9333");
    /// ```
    pub fn weed(master_url: &str, master_port: &str, localhost: bool) -> Self {
        Self {
            master_url: master_url.to_string(),
            master_port: master_port.to_string(),
            localhost,
            ..Default::default()
        }
    }

    /// Build a configuration from environment variables, falling back to
    /// localhost defaults.
    ///
    /// | Variable | Field | Default |
    /// |----------|-------|---------|
    /// | `CEREBRO_FS_URL` | `master_url` | `http://localhost` |
    /// | `CEREBRO_FS_PORT` | `master_port` | `9333` |
    /// | `CEREBRO_FS_FILER_URL` | `filer_url` | `http://localhost:8888` |
    /// | `CEREBRO_FS_ACCESS` | `access` | `weed` |
    /// | `CEREBRO_DANGER_ACCEPT_INVALID_TLS_CERTIFICATE` | `danger_invalid_certificate` | `false` |
    ///
    /// `CEREBRO_FS_ACCESS` accepts `weed` or `filer` (case-insensitive); any
    /// other value falls back to the default and is ignored.
    pub fn from_env() -> Self {
        let access = match env::var("CEREBRO_FS_ACCESS")
            .ok()
            .map(|s| s.to_ascii_lowercase())
            .as_deref()
        {
            Some("filer") => FsAccessMode::Filer,
            Some("weed") => FsAccessMode::Weed,
            _ => FsAccessMode::default(),
        };

        Self {
            master_url: env::var("CEREBRO_FS_URL")
                .unwrap_or_else(|_| DEFAULT_MASTER_URL.to_string()),
            master_port: env::var("CEREBRO_FS_PORT")
                .unwrap_or_else(|_| DEFAULT_MASTER_PORT.to_string()),
            filer_url: env::var("CEREBRO_FS_FILER_URL")
                .unwrap_or_else(|_| DEFAULT_FILER_URL.to_string()),
            localhost: true,
            access,
            danger_invalid_certificate: env::var("CEREBRO_DANGER_ACCEPT_INVALID_TLS_CERTIFICATE")
                .map(|v| v == "true" || v == "1")
                .unwrap_or(false),
        }
    }

    /// Resolve the master HTTP base URL used for cluster health and direct
    /// volume operations.
    ///
    /// In localhost/dev mode the port is appended to the host; in production the
    /// routed `master_url` is used verbatim (the port is assumed to be handled
    /// by the reverse proxy). This preserves the semantics of the previous
    /// `FileSystemClient::get_url` implementation.
    pub fn master_health_url(&self) -> String {
        if self.localhost {
            format!("{}:{}", self.master_url, self.master_port)
        } else {
            self.master_url.clone()
        }
    }

    /// Filer base URL with any trailing slash removed.
    pub fn filer_base(&self) -> String {
        self.filer_url.trim_end_matches('/').to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn master_health_url_localhost_appends_port() {
        let cfg = FsConfig::weed("http://localhost", "9333", true);
        assert_eq!(cfg.master_health_url(), "http://localhost:9333");
    }

    #[test]
    fn master_health_url_routed_uses_url_verbatim() {
        let cfg = FsConfig::weed("https://fs.example.org", "9333", false);
        assert_eq!(cfg.master_health_url(), "https://fs.example.org");
    }

    #[test]
    fn filer_base_trims_trailing_slash() {
        let mut cfg = FsConfig::default();
        cfg.filer_url = "http://localhost:8888/".to_string();
        assert_eq!(cfg.filer_base(), "http://localhost:8888");
    }

    #[test]
    fn access_mode_display_matches_cli_values() {
        assert_eq!(FsAccessMode::Weed.to_string(), "weed");
        assert_eq!(FsAccessMode::Filer.to_string(), "filer");
    }

    #[test]
    fn default_access_is_weed() {
        assert_eq!(FsConfig::default().access, FsAccessMode::Weed);
    }
}
