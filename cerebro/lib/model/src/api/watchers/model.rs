use serde::{Serialize, Deserialize};

use super::schema::RegisterWatcherSchema;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProductionWatcher {
    pub id: String,
    pub date: String,
    pub name: String,
    pub location: String,
    pub last_ping: String,
}
impl ProductionWatcher {
    pub fn from_schema(schema: &RegisterWatcherSchema) -> Self {
        Self {
            id: schema.id.clone(),
            date: schema.date.clone(),
            name: schema.name.clone(),
            location: schema.location.clone(),
            last_ping: schema.last_ping.clone()
        }
    }
}

