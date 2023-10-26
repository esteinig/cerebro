
use serde::{Deserialize, Serialize};

use crate::api::{users::model::UserId, config::Config};

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RegisterTeamSchema {
    pub team_name: String,
    pub team_description: String,
    pub team_lead: UserId,
    pub database_name: String,
    pub database_description: String,
    pub database_mongo_name: String,
    pub project_name: String,
    pub project_description: String,
    pub project_mongo_name: String
}
impl RegisterTeamSchema {
    // A high level methods to verify that 
    // the collection name does not correspond
    // to the logs collection for this team 
    pub fn is_valid(&self, app_config: &Config) -> bool {
        self.project_mongo_name != app_config.database.names.team_database_logs_collection &&
        self.project_mongo_name != app_config.database.names.team_database_reports_collection
    }
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpdateTeamSchema {
    pub team_name: Option<String>,
    pub team_description: String
}

#[derive(Debug, Deserialize, Serialize)]
pub struct RegisterDatabaseSchema {
    pub database_name: String,
    pub database_description: String,
    pub database_mongo_name: String,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct RegisterProjectSchema {
    pub project_name: String,
    pub project_description: String,
    pub project_mongo_name: String,
}