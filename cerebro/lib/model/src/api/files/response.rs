use crate::api::files::model::SeaweedFile;
use serde::{Deserialize, Serialize};

use super::model::FileTag;

#[derive(Serialize, Deserialize)]
pub struct RegisterFileResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>,
}
impl RegisterFileResponse {
    pub fn success(id: &str) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File registered successfully"),
            data: Some(id.to_string()),
        }
    }
    pub fn conflict(id: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("File with unique identifier already exists in database"),
            data: Some(id.to_string()),
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database operation: {}", error_message),
            data: None,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct ListFilesResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<SeaweedFile>>,
}
impl ListFilesResponse {
    pub fn success(files: Vec<SeaweedFile>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File entries found in database"),
            data: Some(files),
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No file entries found in database"),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct UpdateTagsResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<FileTag>>,
}
impl UpdateTagsResponse {
    pub fn success() -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Tags updated for files"),
            data: None,
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Tags updated for files"),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

/// Result data for an expiry (quarantine) sweep.
#[derive(Serialize, Deserialize)]
pub struct ExpiryData {
    /// Number of files quarantined (or eligible, when `dry_run`).
    pub quarantined: u64,
    pub dry_run: bool,
}

/// Response for the non-destructive expiry sweep (`POST /files/expire`).
#[derive(Serialize, Deserialize)]
pub struct ExpireResponse {
    pub status: String,
    pub message: String,
    pub data: Option<ExpiryData>,
}
impl ExpireResponse {
    pub fn success(quarantined: u64, dry_run: bool) -> Self {
        let verb = if dry_run {
            "eligible for quarantine"
        } else {
            "quarantined"
        };
        Self {
            status: String::from("success"),
            message: format!("{} file(s) {}", quarantined, verb),
            data: Some(ExpiryData {
                quarantined,
                dry_run,
            }),
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

/// Result data for a hard purge.
#[derive(Serialize, Deserialize)]
pub struct PurgeData {
    /// Number of files purged (or eligible, when `dry_run`).
    pub purged: u64,
    /// Storage `fid`s of the purged files, for downstream byte reclamation.
    pub fids: Vec<String>,
    pub dry_run: bool,
}

/// Response for the gated hard purge (`POST /files/purge`).
#[derive(Serialize, Deserialize)]
pub struct PurgeResponse {
    pub status: String,
    pub message: String,
    pub data: Option<PurgeData>,
}
impl PurgeResponse {
    pub fn success(fids: Vec<String>, dry_run: bool) -> Self {
        let purged = fids.len() as u64;
        let verb = if dry_run {
            "eligible for purge"
        } else {
            "purged"
        };
        Self {
            status: String::from("success"),
            message: format!("{} file(s) {}", purged, verb),
            data: Some(PurgeData {
                purged,
                fids,
                dry_run,
            }),
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

/// Data for a restore transition: the resulting state and (when restored) the
/// availability-window expiry.
#[derive(Serialize, Deserialize)]
pub struct RestoreTransitionData {
    pub restore_state: crate::api::files::retention::RestoreState,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub restore_expires_at: Option<chrono::DateTime<chrono::Utc>>,
}

/// Response for the restore transition endpoint (`POST /files/{id}/restore`).
#[derive(Serialize, Deserialize)]
pub struct RestoreTransitionResponse {
    pub status: String,
    pub message: String,
    pub data: Option<RestoreTransitionData>,
}
impl RestoreTransitionResponse {
    pub fn success(
        restore_state: crate::api::files::retention::RestoreState,
        restore_expires_at: Option<chrono::DateTime<chrono::Utc>>,
    ) -> Self {
        Self {
            status: String::from("success"),
            message: format!("Restore state is now {}", restore_state),
            data: Some(RestoreTransitionData {
                restore_state,
                restore_expires_at,
            }),
        }
    }
    pub fn noop(restore_state: crate::api::files::retention::RestoreState) -> Self {
        Self {
            status: String::from("success"),
            message: format!("Restore state already {}; no change", restore_state),
            data: Some(RestoreTransitionData {
                restore_state,
                restore_expires_at: None,
            }),
        }
    }
    pub fn invalid_transition(from: &str, to: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("Invalid restore transition {} -> {}", from, to),
            data: None,
        }
    }
    pub fn not_found(id: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("No file found for id {}", id),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

/// Data for an audit-trail query: the matching events and whether the full chain
/// verified intact.
#[derive(Serialize, Deserialize)]
pub struct AuditTrailData {
    pub events: Vec<crate::api::files::audit::AuditEvent>,
    /// Result of verifying the entire team chain's integrity.
    pub verified: bool,
}

/// Response for the audit-trail endpoint (`GET /audit`).
#[derive(Serialize, Deserialize)]
pub struct AuditTrailResponse {
    pub status: String,
    pub message: String,
    pub data: Option<AuditTrailData>,
}
impl AuditTrailResponse {
    pub fn success(events: Vec<crate::api::files::audit::AuditEvent>, verified: bool) -> Self {
        Self {
            status: String::from("success"),
            message: format!("Audit trail retrieved ({} event(s))", events.len()),
            data: Some(AuditTrailData { events, verified }),
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No audit events found"),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

/// Response for fetching a single file by id (`GET /files/{id}`).
#[derive(Serialize, Deserialize)]
pub struct GetFileResponse {
    pub status: String,
    pub message: String,
    pub data: Option<SeaweedFile>,
}
impl GetFileResponse {
    pub fn success(file: SeaweedFile) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File retrieved"),
            data: Some(file),
        }
    }
    pub fn not_found(id: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("No file found for id {}", id),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

/// Result data for a report-out trigger.
#[derive(Serialize, Deserialize)]
pub struct ReportOutData {
    /// Number of files transitioned.
    pub updated: u64,
    /// The report-out timestamp applied (RFC 3339).
    pub reported_at: String,
}

/// Response for the report-out trigger (`POST /files/report-out`).
#[derive(Serialize, Deserialize)]
pub struct ReportOutResponse {
    pub status: String,
    pub message: String,
    pub data: Option<ReportOutData>,
}
impl ReportOutResponse {
    pub fn success(updated: u64, reported_at: String) -> Self {
        Self {
            status: String::from("success"),
            message: format!("Reported out {} file(s)", updated),
            data: Some(ReportOutData {
                updated,
                reported_at,
            }),
        }
    }
    pub fn not_found(run_id: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("No files found for run {}", run_id),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

/// Response for the lifecycle update endpoint (`PATCH /files/{id}/lifecycle`).
#[derive(Serialize, Deserialize)]
pub struct UpdateLifecycleResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>,
}
impl UpdateLifecycleResponse {
    pub fn success(id: &str) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File lifecycle updated"),
            data: Some(id.to_string()),
        }
    }
    pub fn not_found(id: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("No file found for id {}", id),
            data: None,
        }
    }
    pub fn invalid_query() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Lifecycle update contained no fields to change"),
            data: None,
        }
    }
    /// The `expected_tier` compare-and-set precondition did not match the file's
    /// current tier (S2-10) — the update was not applied. A mover treats this as
    /// an idempotent no-op (the move was already handled).
    pub fn precondition_failed(id: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("Tier precondition not met for id {}; update skipped", id),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

/// Represents the output of a successful `weed` command execution.
#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WeedUploadResponse {
    pub file_name: String,
    pub url: String,
    pub fid: String,
    pub size: u64,
}

#[derive(Serialize, Deserialize)]
pub struct DeleteFileResponse {
    pub status: String,
    pub message: String,
    pub data: Option<SeaweedFile>,
}
impl DeleteFileResponse {
    pub fn success(file: SeaweedFile) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File entries deleted from database"),
            data: Some(file),
        }
    }
    /// Deletion refused because the file is under legal hold or still within its
    /// retention period (S2-9). Clear the hold / retention via the lifecycle
    /// endpoint first (an audited admin action).
    pub fn protected(reason: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("File is protected from deletion: {}", reason),
            data: None,
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No file entries found in database"),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct DeleteFilesResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<String>>,
}
impl DeleteFilesResponse {
    pub fn success(fids: Vec<String>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("File entries deleted from database"),
            data: Some(fids),
        }
    }
    pub fn invalid_query() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Invalid query - run identifier or sample identifier must be provided if no query parameter 'all' is porovided"),
            data: None
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No file entries found in database"),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None,
        }
    }
}
