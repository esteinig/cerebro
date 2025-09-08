use serde::{Deserialize, Serialize};
use crate::api::cerebro::schema::PrefetchData;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateTrainingPrefetch {
    pub id: String,
    pub collection: String,
    pub identifier: String,
    pub name: String,
    pub prefetch: PrefetchData,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryTrainingPrefetch {
    /// Optional filter by collection
    pub collection: Option<String>,
    /// Optional filter by identifier
    pub identifier: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeleteTrainingPrefetch {
    pub id: String, // hex ObjectId
}