use serde::{Deserialize, Serialize};

use crate::api::training::schema::CreateTrainingPrefetch;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingPrefetchRecord {
    // Unique identifier for retrieving prefetch data from GridFs
    pub id: String,
    /// Name of a training collection
    pub collection: String,
    /// Unique identifier within the collection
    pub identifier: String,
    /// Human readable label
    pub name: String,
}
impl TrainingPrefetchRecord {
    pub fn from_request(req: &CreateTrainingPrefetch) -> Self {
        Self {
            id: req.id.clone(),
            collection: req.collection.clone(),
            identifier: req.identifier.clone(),
            name: req.name.clone()
        }
    }
}