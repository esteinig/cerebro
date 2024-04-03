use serde::{Serialize, Deserialize};
use crate::api::files::model::{SeaweedFileId, WatcherConfig};

#[derive(Deserialize, Serialize)]
pub struct RegisterFileSchema {
    pub id: String,
    pub date: String,
    pub name: String,
    pub hash: String,
    pub fid: SeaweedFileId,
    pub size: u64,
    pub watcher: WatcherConfig
}
