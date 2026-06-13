//! Integrity verification and replica repair for Cerebro FS.
//!
//! Files are registered with a BLAKE3 hash at upload. This module sweeps
//! registered files, re-reads their bytes, and confirms the hash still matches —
//! the core of an accreditation-grade integrity check. On a mismatch it can
//! **repair** the local copy by re-fetching from an alternate SeaweedFS replica
//! (a fid is served by every volume server that holds a copy), and it reports
//! the observed replica count per file.
//!
//! Archived (Glacier) objects are skipped: they must be restored first.
//!
//! ## SeaweedFS endpoints
//!
//! Replica discovery uses the master volume lookup
//! (`GET /dir/lookup?volumeId=<id>`), and a per-replica fetch is a direct
//! `GET http://<volume>/<fid>`. Confirm these against your SeaweedFS version
//! (the stack pins 3.64).

use std::path::PathBuf;
use std::time::Duration;

use reqwest::StatusCode;
use serde::Deserialize;

use cerebro_model::api::files::model::SeaweedFile;

use crate::client::{verify_file, FileSystemClient};
use crate::error::FileSystemError;
use crate::hash::fast_file_hash;

/// Result of verifying a single file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VerifyStatus {
    /// The recomputed hash matched the registered hash.
    Verified,
    /// The first copy mismatched but an alternate replica verified successfully.
    Repaired,
    /// The hash did not match and no good replica was found.
    Mismatch { expected: String, actual: String },
    /// The file could not be retrieved.
    Missing,
    /// Verification was not attempted (e.g. the object is archived).
    Skipped(String),
}

/// Verification outcome for one file, including the observed replica count.
#[derive(Debug, Clone)]
pub struct VerifyOutcome {
    pub identifier: String,
    pub name: String,
    pub replicas: Option<u32>,
    pub status: VerifyStatus,
}

/// Aggregate report of an integrity sweep.
#[derive(Debug, Clone)]
pub struct VerifyReport {
    pub outcomes: Vec<VerifyOutcome>,
}
impl VerifyReport {
    /// Number of files that verified cleanly or were repaired.
    pub fn ok_count(&self) -> usize {
        self.outcomes
            .iter()
            .filter(|o| matches!(o.status, VerifyStatus::Verified | VerifyStatus::Repaired))
            .count()
    }
    /// Number of files that failed verification (mismatch or missing).
    pub fn failed_count(&self) -> usize {
        self.outcomes
            .iter()
            .filter(|o| matches!(o.status, VerifyStatus::Mismatch { .. } | VerifyStatus::Missing))
            .count()
    }
    /// True when no file failed verification.
    pub fn ok(&self) -> bool {
        self.failed_count() == 0
    }
}

/// A SeaweedFS volume location returned by the master lookup.
#[derive(Debug, Clone, Deserialize)]
pub struct VolumeLocation {
    pub url: String,
    #[serde(default, rename = "publicUrl")]
    pub public_url: String,
}

#[derive(Debug, Clone, Deserialize)]
struct LookupResponse {
    #[serde(default)]
    locations: Vec<VolumeLocation>,
}

/// Extract the SeaweedFS volume id from a fid of the form `<volumeId>,<fileKey>`.
pub fn volume_id_from_fid(fid: &str) -> Option<&str> {
    match fid.split_once(',') {
        Some((volume_id, _)) if !volume_id.is_empty() => Some(volume_id),
        _ => None,
    }
}

impl FileSystemClient {
    /// Look up the volume locations (replicas) that serve a given fid.
    ///
    /// Filer paths have no single fid and return an empty list.
    pub fn lookup_locations(&self, fid: &str) -> Result<Vec<VolumeLocation>, FileSystemError> {
        let volume_id = match volume_id_from_fid(fid) {
            Some(id) => id,
            None => return Ok(Vec::new()),
        };

        let url = format!("{}/dir/lookup?volumeId={}", self.get_url(), volume_id);
        let response = reqwest::blocking::Client::new().get(&url).send()?;

        match response.status() {
            StatusCode::OK => {
                let parsed: LookupResponse = response.json()?;
                Ok(parsed.locations)
            }
            status => Err(FileSystemError::UnexpectedResponseStatus(status)),
        }
    }

    /// Download a fid directly from a specific volume server (used for repair).
    fn download_from_volume(
        &self,
        location: &VolumeLocation,
        fid: &str,
        out_path: &PathBuf,
    ) -> Result<(), FileSystemError> {
        // Volume URLs from the master lack a scheme; assume http on the cluster network.
        let base = if location.url.starts_with("http") {
            location.url.clone()
        } else {
            format!("http://{}", location.url)
        };
        let url = format!("{}/{}", base, fid);

        let http = reqwest::blocking::Client::builder()
            .danger_accept_invalid_certs(self.config.danger_invalid_certificate)
            .timeout(Duration::from_secs(600))
            .build()?;
        let mut response = http.get(&url).send()?;
        if !response.status().is_success() {
            return Err(FileSystemError::UnexpectedResponseStatus(response.status()));
        }
        let mut file = std::fs::File::create(out_path)?;
        response.copy_to(&mut file)?;
        Ok(())
    }

    /// Verify the integrity of files in a run/sample against their registered
    /// BLAKE3 hashes, optionally repairing mismatches from an alternate replica.
    ///
    /// Each non-archived file is downloaded to a temporary directory and hashed.
    /// When `repair` is set and a mismatch is found, every other replica reported
    /// by the master is tried until one verifies. Archived (Glacier) objects are
    /// skipped. The temporary copies are removed before returning.
    pub fn verify_files(
        &self,
        run_id: Option<String>,
        sample_id: Option<String>,
        repair: bool,
    ) -> Result<VerifyReport, FileSystemError> {
        let run_id = run_id.ok_or(FileSystemError::InvalidDownloadQuery)?;

        let files = self.api_client.list_files(Some(run_id), None, 0, 100_000, false)?;
        let files: Vec<SeaweedFile> = match &sample_id {
            Some(sid) => files
                .into_iter()
                .filter(|f| f.sample_id.as_deref() == Some(sid.as_str()))
                .collect(),
            None => files,
        };

        let tmp = std::env::temp_dir().join(format!("cerebro-verify-{}", std::process::id()));
        std::fs::create_dir_all(&tmp)?;

        let mut outcomes = Vec::new();
        for file in &files {
            let identifier = file.effective_identifier().to_string();

            if file.requires_restore() {
                outcomes.push(VerifyOutcome {
                    identifier,
                    name: file.name.clone(),
                    replicas: None,
                    status: VerifyStatus::Skipped("archived; restore required".to_string()),
                });
                continue;
            }

            let replicas = self
                .lookup_locations(&file.fid)
                .map(|l| l.len() as u32)
                .ok();

            let status = self.verify_one(file, &tmp, repair)?;
            outcomes.push(VerifyOutcome {
                identifier,
                name: file.name.clone(),
                replicas,
                status,
            });
        }

        // Best-effort cleanup of the temporary working directory.
        let _ = std::fs::remove_dir_all(&tmp);

        Ok(VerifyReport { outcomes })
    }

    /// Verify a single file, attempting replica repair on mismatch when enabled.
    fn verify_one(
        &self,
        file: &SeaweedFile,
        tmp: &PathBuf,
        repair: bool,
    ) -> Result<VerifyStatus, FileSystemError> {
        let identifier = file.effective_identifier().to_string();
        let target = tmp.join(&file.name);

        // Initial fetch via the normal routing (filer path or weed fid).
        if self
            .download_identifier(&identifier, Some(file.name.as_str()), tmp)
            .is_err()
        {
            return Ok(VerifyStatus::Missing);
        }
        let local = if target.exists() { target.clone() } else { tmp.join(&file.name) };

        match verify_file(&local, &file.hash) {
            Ok(()) => Ok(VerifyStatus::Verified),
            Err(FileSystemError::IntegrityMismatch { expected, actual, .. }) => {
                if !repair {
                    return Ok(VerifyStatus::Mismatch { expected, actual });
                }
                // Repair: try each alternate replica until one verifies.
                let locations = self.lookup_locations(&file.fid).unwrap_or_default();
                for location in &locations {
                    let _ = std::fs::remove_file(&local);
                    if self.download_from_volume(location, &file.fid, &local).is_err() {
                        continue;
                    }
                    if let Ok(actual) = fast_file_hash(&local) {
                        if actual == file.hash {
                            return Ok(VerifyStatus::Repaired);
                        }
                    }
                }
                Ok(VerifyStatus::Mismatch { expected, actual })
            }
            Err(other) => Err(other),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn volume_id_parsed_from_fid() {
        assert_eq!(volume_id_from_fid("3,01637037d6"), Some("3"));
        assert_eq!(volume_id_from_fid("17,abc"), Some("17"));
        assert_eq!(volume_id_from_fid("no-comma"), None);
        assert_eq!(volume_id_from_fid(",01"), None);
    }

    #[test]
    fn lookup_response_parses_locations() {
        let body = r#"{"volumeOrFileId":"3","locations":[{"url":"10.0.0.1:8080","publicUrl":"node1:8080"},{"url":"10.0.0.2:8080","publicUrl":"node2:8080"}]}"#;
        let parsed: LookupResponse = serde_json::from_str(body).unwrap();
        assert_eq!(parsed.locations.len(), 2);
        assert_eq!(parsed.locations[0].url, "10.0.0.1:8080");
        assert_eq!(parsed.locations[1].public_url, "node2:8080");
    }

    #[test]
    fn report_counts_and_ok() {
        let report = VerifyReport {
            outcomes: vec![
                VerifyOutcome { identifier: "a".into(), name: "a".into(), replicas: Some(2), status: VerifyStatus::Verified },
                VerifyOutcome { identifier: "b".into(), name: "b".into(), replicas: Some(2), status: VerifyStatus::Repaired },
                VerifyOutcome { identifier: "c".into(), name: "c".into(), replicas: Some(1), status: VerifyStatus::Mismatch { expected: "x".into(), actual: "y".into() } },
                VerifyOutcome { identifier: "d".into(), name: "d".into(), replicas: None, status: VerifyStatus::Skipped("archived".into()) },
            ],
        };
        assert_eq!(report.ok_count(), 2);
        assert_eq!(report.failed_count(), 1);
        assert!(!report.ok());
    }
}
