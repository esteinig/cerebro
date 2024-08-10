use std::{fs::File, path::PathBuf, io::{Write, Read}};
use serde::{Serialize, Deserialize};
use crate::api::cerebro::model::ModelError;

use super::model::WatcherFormat;

#[derive(Deserialize, Serialize, Debug)]
pub struct RegisterWatcherSchema {
    pub id: String,
    pub date: String,
    pub name: String,
    pub location: String,
    pub format: WatcherFormat,
    pub glob: String,
    pub last_ping: String
}
impl RegisterWatcherSchema {
    pub fn new(name: &str, location: &str, format: WatcherFormat, glob: Option<String>) -> Self {
        
        let date = chrono::Utc::now().to_string();
        let glob = glob.unwrap_or(format.default_glob());

        Self {
            id: uuid::Uuid::new_v4().to_string(),
            date: date.clone(),
            name: name.to_string(),
            location: location.to_string(),
            format,
            glob,
            last_ping: date
        }
    }
    pub fn from_json(path: &PathBuf) -> Result<Self, ModelError> {
        let mut file = File::open(path)?;
        let mut json_data = String::new();
        file.read_to_string(&mut json_data)?;
        let schema = serde_json::from_str(&json_data).map_err(ModelError::JsonDeserialization)?;
        Ok(schema)
    }
    pub fn to_json(&self, path: &PathBuf) -> Result<(), ModelError> {
        let json_data = serde_json::to_string_pretty(&self).map_err(ModelError::JsonSerialization)?;
        let mut file = File::create(path)?;
        file.write_all(json_data.as_bytes())?;

        Ok(())
    }
}