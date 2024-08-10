use serde::{Serialize, Deserialize};
use crate::api::{files::model::SeaweedFileId, pipelines::model::ProductionPipeline, watchers::model::ProductionWatcher};

#[derive(Deserialize, Serialize, Debug)]
pub struct RegisterFileSchema {
    pub id: String,
    pub run_id: Option<String>,
    pub sample_id: Option<String>,
    pub date: String,
    pub name: String,
    pub hash: String,
    pub fid: SeaweedFileId,
    pub size: u64,
    pub watcher: Option<ProductionWatcher>,
    pub pipeline: Option<ProductionPipeline>
}
