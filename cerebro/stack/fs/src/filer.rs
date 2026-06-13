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
//! Directory listing and tiering operations are added in later components.
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
            return Err(FilerError::LocalFileMissing(local_path.display().to_string()));
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
}
