use std::path::Path;

use serde::{Deserialize, Serialize};
use uuid::Uuid;
use crate::api::cerebro::{model::ModelError, schema::PrefetchData};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateTrainingPrefetch {
    pub id: String,
    pub collection: String,
    pub identifier: String,
    pub name: String,
    pub prefetch: PrefetchData,
}
impl CreateTrainingPrefetch {
    pub fn from_file(path: &Path, collection: &str) -> Result<Self, ModelError> {

        let prefetch_data = PrefetchData::from_json(path)?;

        Ok(Self {
            id: Uuid::new_v4().to_string(),
            collection: collection.to_string(),
            identifier: prefetch_data.config.sample.clone(),
            name: prefetch_data.config.sample.clone(),
            prefetch: prefetch_data
        })

    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryTrainingPrefetch {
    /// Optional filter by collection
    pub collection: Option<String>,
    /// Optional filter by identifier
    pub identifier: Option<String>,
    /// Optional filter by name
    pub name: Option<String>
}
