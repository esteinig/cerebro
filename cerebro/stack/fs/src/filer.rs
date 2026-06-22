//! Minimal SeaweedFS Filer HTTP client.
//!
//! The filer exposes a path-addressed HTTP API on top of the SeaweedFS volume
//! servers. Compared with the `weed upload`/`weed download` shell-outs in
//! [`crate::weed`], the filer:
//!
//! * gives objects a stable **path** (e.g. `/run/sample/reads_R1.fastq.gz`)
//!   rather than an opaque fid;
//! * records that path in the filer metadata store (MongoDB in the Cerebro
//!   deployment), enabling directory listings and lifecycle operations; and
//! * auto-chunks large objects server-side, so uploads/downloads can be
//!   **streamed** from/to disk without buffering whole files in memory.
//!
//! Directory listing ([`FilerClient::list_objects`]) backs store-side orphan
//! detection; an object-presence probe ([`FilerClient::exists`]) backs the archival
//! reclaim's safety gate.
//!
//! ## Large files
//!
//! Sequencing read sets are routinely multi-gigabyte. [`FilerClient::upload`]
//! sends the file as a multipart part that is read lazily from disk, and
//! [`FilerClient::download`] streams the response straight to a file, so peak
//! memory use is independent of object size. No global request timeout is set
//! on the HTTP client for this reason.

use std::fs::{create_dir_all, File};
use std::path::Path;
use std::time::Duration;

use chrono::{DateTime, Utc};
use reqwest::blocking::{multipart, Client};
use reqwest::StatusCode;
use serde::Deserialize;

use crate::error::FilerError;

/// Response body returned by the filer for a successful upload.
///
/// All fields are optional: for small objects the filer returns a single `fid`,
/// while large (auto-chunked) objects are addressed purely by their path and may
/// omit it. Callers should therefore treat the upload **path** as the durable
/// identifier and use `size` only as an advisory value.
#[derive(Debug, Clone, Deserialize)]
pub struct FilerUploadResponse {
    /// Stored file name, if reported.
    pub name: Option<String>,
    /// Object size in bytes, if reported.
    pub size: Option<u64>,
    /// Assigned fid for non-chunked objects, if reported.
    pub fid: Option<String>,
    /// Error message reported by the filer, if any.
    pub error: Option<String>,
}

/// HTTP client for a single SeaweedFS filer.
#[derive(Clone, Debug)]
pub struct FilerClient {
    base_url: String,
    http: Client,
}

impl FilerClient {
    /// Create a new filer client for the given base URL
    /// (e.g. `http://localhost:8888`).
    ///
    /// No request timeout is configured: large sequencing uploads can take many
    /// minutes and a global timeout would abort them. A connect timeout still
    /// guards against an unreachable host.
    pub fn new(base_url: &str, danger_invalid_certificate: bool) -> Result<Self, FilerError> {
        let http = Client::builder()
            .danger_accept_invalid_certs(danger_invalid_certificate)
            .connect_timeout(Duration::from_secs(30))
            .build()?;

        Ok(Self {
            base_url: base_url.trim_end_matches('/').to_string(),
            http,
        })
    }

    /// Build the absolute object URL for a remote path, normalising slashes so
    /// that exactly one separator joins the base URL and the path.
    fn object_url(&self, remote_path: &str) -> String {
        format!("{}/{}", self.base_url, remote_path.trim_start_matches('/'))
    }

    /// Upload a local file to `remote_path`, streaming the body from disk.
    ///
    /// The file is sent as a multipart form part read lazily from disk, so peak
    /// memory use is independent of file size. The filer assigns volume chunks
    /// as needed for large objects.
    ///
    /// # Errors
    ///
    /// Returns [`FilerError::LocalFileMissing`] if `local_path` does not exist,
    /// [`FilerError::UnexpectedStatus`] for a non-success HTTP status, or
    /// [`FilerError::Upload`] if the filer reports an application-level error.
    pub fn upload(
        &self,
        local_path: &Path,
        remote_path: &str,
    ) -> Result<FilerUploadResponse, FilerError> {
        if !local_path.exists() {
            return Err(FilerError::LocalFileMissing(
                local_path.display().to_string(),
            ));
        }

        // `multipart::Form::file` opens the file and streams it; the whole file
        // is never resident in memory.
        let form = multipart::Form::new().file("file", local_path)?;

        let url = self.object_url(remote_path);
        let response = self.http.post(&url).multipart(form).send()?;

        let status = response.status();
        if !status.is_success() {
            return Err(FilerError::UnexpectedStatus(status));
        }

        // The filer returns JSON on success; tolerate an empty body by
        // synthesising an (otherwise empty) response.
        let body = response.text()?;
        if body.trim().is_empty() {
            return Ok(FilerUploadResponse {
                name: None,
                size: None,
                fid: None,
                error: None,
            });
        }

        let parsed: FilerUploadResponse = serde_json::from_str(&body)?;
        if let Some(err) = &parsed.error {
            if !err.is_empty() {
                return Err(FilerError::Upload(err.clone()));
            }
        }
        Ok(parsed)
    }

    /// Download `remote_path` to `out_path`, streaming the response to disk.
    ///
    /// Parent directories of `out_path` are created if necessary.
    ///
    /// # Errors
    ///
    /// Returns [`FilerError::NotFound`] if the object does not exist, or
    /// [`FilerError::UnexpectedStatus`] for any other non-`200` status.
    pub fn download(&self, remote_path: &str, out_path: &Path) -> Result<(), FilerError> {
        if let Some(parent) = out_path.parent() {
            if !parent.as_os_str().is_empty() && !parent.exists() {
                create_dir_all(parent)?;
            }
        }

        let url = self.object_url(remote_path);
        let mut response = self.http.get(&url).send()?;

        match response.status() {
            StatusCode::OK => {
                let mut file = File::create(out_path)?;
                // `copy_to` streams the body in chunks rather than buffering it.
                response.copy_to(&mut file)?;
                Ok(())
            }
            StatusCode::NOT_FOUND => Err(FilerError::NotFound(remote_path.to_string())),
            status => Err(FilerError::UnexpectedStatus(status)),
        }
    }

    /// Stream `remote_path` from the filer and return its BLAKE3 digest, without
    /// writing it to disk.
    ///
    /// The response body implements [`Read`](std::io::Read), so it is piped
    /// directly into a `blake3::Hasher` in constant memory — suitable for the
    /// multi-gigabyte read sets this platform stores. No temp file is created.
    pub fn hash(&self, remote_path: &str) -> Result<String, FilerError> {
        let url = self.object_url(remote_path);
        let mut response = self.http.get(&url).send()?;

        match response.status() {
            StatusCode::OK => {
                let mut hasher = blake3::Hasher::new();
                // Streams the body in chunks rather than buffering the whole object.
                std::io::copy(&mut response, &mut hasher)?;
                Ok(hasher.finalize().to_hex().to_string())
            }
            StatusCode::NOT_FOUND => Err(FilerError::NotFound(remote_path.to_string())),
            status => Err(FilerError::UnexpectedStatus(status)),
        }
    }

    /// Delete `remote_path` from the filer.
    ///
    /// When `recursive` is set, a directory and all of its contents are removed.
    pub fn delete(&self, remote_path: &str, recursive: bool) -> Result<(), FilerError> {
        let mut url = self.object_url(remote_path);
        if recursive {
            url.push_str("?recursive=true");
        }

        let response = self.http.delete(&url).send()?;
        match response.status() {
            StatusCode::OK | StatusCode::ACCEPTED | StatusCode::NO_CONTENT => Ok(()),
            StatusCode::NOT_FOUND => Err(FilerError::NotFound(remote_path.to_string())),
            status => Err(FilerError::UnexpectedStatus(status)),
        }
    }

    /// Lightweight reachability check: a successful GET of the filer root is
    /// treated as healthy.
    ///
    /// A full topology health check (master + filer + volumes + S3); this method
    /// exists so callers can fail fast on an unreachable filer.
    pub fn health(&self) -> Result<(), FilerError> {
        let url = format!("{}/", self.base_url);
        let response = self.http.get(&url).send()?;
        if response.status().is_success() {
            Ok(())
        } else {
            Err(FilerError::UnexpectedStatus(response.status()))
        }
    }

    /// Whether an object exists at `remote_path`: a cheap `HEAD`. `200` =>
    /// present, `404` => absent. Used by the archival reclaim to confirm a local
    /// filer copy is present before deleting it.
    pub fn exists(&self, remote_path: &str) -> Result<bool, FilerError> {
        let url = self.object_url(remote_path);
        let response = self.http.head(&url).send()?;
        match response.status() {
            StatusCode::NOT_FOUND => Ok(false),
            s if s.is_success() => Ok(true),
            status => Err(FilerError::UnexpectedStatus(status)),
        }
    }

    /// Recursively list **file** object paths under `base`, bounded by
    /// `max_objects`.
    ///
    /// Walks the filer directory tree depth-first, paginating each directory, and
    /// returns one [`FilerObject`] per file (directories are traversed, not
    /// returned). The total object budget caps the walk so a pathologically large
    /// estate can never run unbounded; when the budget is hit the partial result is
    /// returned (the caller treats a truncated listing as inconclusive for orphan
    /// deletion). Used by the consistency-reconcile scan to enumerate what the
    /// store actually holds, for orphan detection.
    pub fn list_objects(
        &self,
        base: &str,
        max_objects: usize,
    ) -> Result<Vec<FilerObject>, FilerError> {
        let mut out: Vec<FilerObject> = Vec::new();
        let mut stack: Vec<String> = vec![base.to_string()];

        while let Some(dir) = stack.pop() {
            if out.len() >= max_objects {
                break;
            }
            let mut last: Option<String> = None;
            loop {
                let page = self.list_dir_page(&dir, last.as_deref())?;
                let entries = page.entries.unwrap_or_default();
                if entries.is_empty() {
                    break;
                }
                let mut last_name = None;
                for e in &entries {
                    last_name = e.full_path.rsplit('/').next().map(|s| s.to_string());
                    if e.is_dir() {
                        stack.push(e.full_path.clone());
                    } else {
                        out.push(FilerObject {
                            path: e.full_path.clone(),
                            mtime: e.parsed_mtime(),
                        });
                        if out.len() >= max_objects {
                            return Ok(out);
                        }
                    }
                }
                if page.should_display_load_more == Some(true) {
                    // Continue the same directory from the last seen name.
                    last = page.last_file_name.clone().or(last_name);
                    if last.is_none() {
                        break;
                    }
                } else {
                    break;
                }
            }
        }
        Ok(out)
    }

    /// Fetch one page of a directory listing as JSON.
    fn list_dir_page(&self, dir: &str, last: Option<&str>) -> Result<FilerListing, FilerError> {
        let mut url = format!("{}/{}", self.base_url, dir.trim_start_matches('/'));
        if !url.ends_with('/') {
            url.push('/');
        }
        url.push_str(&format!("?limit={LIST_PAGE_LIMIT}"));
        if let Some(l) = last {
            url.push_str("&lastFileName=");
            url.push_str(&l.replace(' ', "%20"));
        }

        let response = self
            .http
            .get(&url)
            .header("Accept", "application/json")
            .send()?;
        match response.status() {
            s if s.is_success() => {
                let body = response.text()?;
                Ok(serde_json::from_str::<FilerListing>(&body)?)
            }
            StatusCode::NOT_FOUND => Ok(FilerListing::default()),
            status => Err(FilerError::UnexpectedStatus(status)),
        }
    }
}

/// Per-directory page size for filer listings.
const LIST_PAGE_LIMIT: u32 = 1000;

/// A file object discovered by a filer listing. In filer mode the **path** is
/// the object's stable identifier (compared against the catalogue's `path`).
#[derive(Debug, Clone)]
pub struct FilerObject {
    pub path: String,
    /// Last-modified time, when the filer reported a parseable RFC3339 timestamp.
    pub mtime: Option<DateTime<Utc>>,
}

/// One entry in a filer directory listing.
#[derive(Debug, Clone, Deserialize)]
struct RawFilerEntry {
    #[serde(rename = "FullPath")]
    full_path: String,
    #[serde(rename = "Mtime")]
    mtime: Option<String>,
    /// Go `os.FileMode`; the directory bit is `1 << 31`.
    #[serde(rename = "Mode", default)]
    mode: u64,
}
impl RawFilerEntry {
    fn is_dir(&self) -> bool {
        self.mode & (1 << 31) != 0
    }
    fn parsed_mtime(&self) -> Option<DateTime<Utc>> {
        self.mtime
            .as_deref()
            .and_then(|s| DateTime::parse_from_rfc3339(s).ok())
            .map(|dt| dt.with_timezone(&Utc))
    }
}

/// A filer directory-listing response page.
#[derive(Debug, Default, Deserialize)]
struct FilerListing {
    #[serde(rename = "Entries")]
    entries: Option<Vec<RawFilerEntry>>,
    #[serde(rename = "LastFileName")]
    last_file_name: Option<String>,
    #[serde(rename = "ShouldDisplayLoadMore")]
    should_display_load_more: Option<bool>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn client() -> FilerClient {
        // Trailing slash on the base URL is intentional: it must be trimmed.
        FilerClient::new("http://localhost:8888/", false).expect("client builds")
    }

    #[test]
    fn object_url_joins_with_single_slash() {
        let c = client();
        assert_eq!(
            c.object_url("run/sample/file.gz"),
            "http://localhost:8888/run/sample/file.gz"
        );
    }

    #[test]
    fn object_url_strips_leading_slash_on_path() {
        let c = client();
        assert_eq!(
            c.object_url("/run/sample/file.gz"),
            "http://localhost:8888/run/sample/file.gz"
        );
    }

    #[test]
    fn object_url_empty_path_yields_base_with_slash() {
        let c = client();
        assert_eq!(c.object_url(""), "http://localhost:8888/");
    }

    // --- Filer listing parser ---

    #[test]
    fn entry_directory_bit_detected() {
        // The directory bit is Go's os.ModeDir = 1 << 31.
        let dir = RawFilerEntry {
            full_path: "/a/b".into(),
            mtime: None,
            mode: 1 << 31,
        };
        assert!(dir.is_dir());
        let dir_with_perms = RawFilerEntry {
            full_path: "/a/b".into(),
            mtime: None,
            mode: (1 << 31) | 0o755,
        };
        assert!(dir_with_perms.is_dir());
        let file = RawFilerEntry {
            full_path: "/a/b.gz".into(),
            mtime: None,
            mode: 0o644,
        };
        assert!(!file.is_dir());
    }

    #[test]
    fn entry_mtime_parses_rfc3339_else_none() {
        let ok = RawFilerEntry {
            full_path: "/f".into(),
            mtime: Some("2024-01-02T03:04:05Z".into()),
            mode: 0,
        };
        assert!(ok.parsed_mtime().is_some());
        let bad = RawFilerEntry {
            full_path: "/f".into(),
            mtime: Some("not-a-date".into()),
            mode: 0,
        };
        assert!(bad.parsed_mtime().is_none()); // unparseable -> None, never an error
        let missing = RawFilerEntry {
            full_path: "/f".into(),
            mtime: None,
            mode: 0,
        };
        assert!(missing.parsed_mtime().is_none());
    }

    #[test]
    fn listing_deserializes_entries_and_pagination() {
        let json = r#"{
            "Entries": [
                { "FullPath": "/run/sub", "Mtime": "2024-01-02T03:04:05Z", "Mode": 2147483648 },
                { "FullPath": "/run/reads.gz", "Mtime": "2024-01-02T03:04:05Z", "Mode": 420 }
            ],
            "LastFileName": "reads.gz",
            "ShouldDisplayLoadMore": true
        }"#;
        let parsed: FilerListing = serde_json::from_str(json).expect("valid listing");
        let entries = parsed.entries.expect("entries present");
        assert_eq!(entries.len(), 2);
        assert!(entries[0].is_dir()); // /run/sub (Mode 2147483648 = 1<<31)
        assert!(!entries[1].is_dir()); // /run/reads.gz (Mode 420 = 0o644)
        assert_eq!(parsed.last_file_name.as_deref(), Some("reads.gz"));
        assert_eq!(parsed.should_display_load_more, Some(true));
    }

    #[test]
    fn listing_tolerates_missing_optionals_and_default_is_empty() {
        // Mode absent -> default 0 (not a dir); optional listing fields absent.
        let json = r#"{ "Entries": [ { "FullPath": "/f" } ] }"#;
        let parsed: FilerListing = serde_json::from_str(json).expect("valid");
        let entries = parsed.entries.expect("entries");
        assert_eq!(entries.len(), 1);
        assert!(!entries[0].is_dir());
        assert!(entries[0].parsed_mtime().is_none());
        assert!(parsed.last_file_name.is_none());
        assert!(parsed.should_display_load_more.is_none());

        // The Default (returned on a 404 page) carries no entries.
        let empty = FilerListing::default();
        assert!(empty.entries.is_none());
    }
}
