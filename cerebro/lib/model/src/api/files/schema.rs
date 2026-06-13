use serde::{Serialize, Deserialize};
use chrono::{DateTime, Utc};
use crate::api::{files::model::SeaweedFileId, watchers::model::ProductionWatcher};

use super::model::{FileTag, FileType};
use super::retention::{RetentionClass, StorageTier};

#[derive(Deserialize, Serialize, Debug)]
pub struct RegisterFileSchema {
    pub id: String,
    pub date: String,
    pub name: String,
    pub hash: String,
    pub size: u64,
    pub fid: SeaweedFileId,
    pub run_id: Option<String>,
    pub sample_id: Option<String>,
    pub pipeline_id: Option<String>,
    pub description: Option<String>,
    pub ftype: Option<FileType>,
    pub watcher: Option<ProductionWatcher>,
    /// Filer object path when stored via the path-addressed filer (FS-2).
    #[serde(default)]
    pub path: Option<String>,
    /// Physical storage tier at registration.
    #[serde(default)]
    pub tier: StorageTier,
    /// Retention category assigned to the file.
    #[serde(default)]
    pub retention: RetentionClass,
    /// Absolute expiry computed by the registering client from its retention
    /// policy. `None` means "retain indefinitely".
    #[serde(default)]
    pub retain_until: Option<DateTime<Utc>>,
    /// When set, the file is exempt from expiry.
    #[serde(default)]
    pub legal_hold: bool,
    /// Requested/observed replica count, when known.
    #[serde(default)]
    pub replicas: Option<u32>,
    /// Whether the object's data resides in remote archival storage (Glacier)
    /// and requires a restore before retrieval.
    #[serde(default)]
    pub archived: bool,
    /// When the result for this file's case was reported out (retention anchor).
    #[serde(default)]
    pub reported_at: Option<DateTime<Utc>>,
}

#[derive(Deserialize, Serialize, Debug)]
pub struct UpdateFileTagsSchema {
    pub ids: Vec<String>,
    pub tags: Vec<FileTag>
}