use serde::{Serialize, Deserialize};

use crate::api::{pipelines::model::ProductionPipeline, watchers::model::ProductionWatcher};

use super::schema::RegisterFileSchema;

/*
========================
File system and storage
========================
*/

pub type SeaweedFileId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeaweedFile {
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
impl SeaweedFile {
    pub fn from_schema(register_file_schema: &RegisterFileSchema) -> Self {
        Self {
            id: register_file_schema.id.clone(),
            run_id: register_file_schema.run_id.clone(),
            sample_id: register_file_schema.sample_id.clone(),
            date: register_file_schema.date.clone(),
            name: register_file_schema.name.clone(),
            hash: register_file_schema.hash.clone(),
            fid: register_file_schema.fid.clone(),
            size: register_file_schema.size.clone(),
            watcher: register_file_schema.watcher.clone(),
            pipeline: register_file_schema.pipeline.clone(),
        }
    }
    pub fn size_mb(&self) -> f64 {
        bytes_to_mb(self.size)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeaweedReads {
    pub reads_1: SeaweedFile,
    pub reads_2: Option<SeaweedFile>
}


/// Converts a size in bytes to megabytes (MB).
///
/// # Arguments
///
/// * `bytes` - The size in bytes.
///
/// # Returns
///
/// The size in megabytes (MB) as a floating point number.
fn bytes_to_mb(bytes: u64) -> f64 {
    const BYTES_PER_MB: f64 = 1024.0 * 1024.0; // 1 MB = 1024 * 1024 bytes
    bytes as f64 / BYTES_PER_MB
}
