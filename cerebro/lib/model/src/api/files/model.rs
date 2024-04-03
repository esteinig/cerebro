use serde::{Serialize, Deserialize};

use super::schema::RegisterFileSchema;

/*
========================
File system and storage
========================
*/

pub type SeaweedFileId = String;


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WatcherConfig {
    pub id: String,
    pub name: String,
    pub location: String
}

impl Default for WatcherConfig {
    fn default() -> Self {
        Self {
            id: uuid::Uuid::new_v4().to_string(),
            name: "Eye of Sauron".to_string(),
            location: "Barad-dûr".to_string()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeaweedFile {
    pub id: String,
    pub date: String,
    pub name: String,
    pub hash: String,
    pub fid: SeaweedFileId,
    pub size: u64,
    pub watcher: WatcherConfig
}
impl SeaweedFile {
    pub fn from_schema(register_file_schema: &RegisterFileSchema) -> Self {
        Self {
            id: register_file_schema.id.clone(),
            date: register_file_schema.date.clone(),
            name: register_file_schema.name.clone(),
            hash: register_file_schema.hash.clone(),
            fid: register_file_schema.fid.clone(),
            size: register_file_schema.size.clone(),
            watcher: register_file_schema.watcher.clone(),
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
