use serde::{Serialize, Deserialize};
use crate::api::{files::model::SeaweedFileId, watchers::model::ProductionWatcher};

use super::model::{FileTag, FileType};

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
    pub ftype: Option<FileType>,
    pub watcher: Option<ProductionWatcher>
}

#[derive(Deserialize, Serialize, Debug)]
pub struct UpdateFileTagsSchema {
    pub ids: Vec<String>,
    pub tags: Vec<FileTag>
}