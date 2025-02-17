
use serde::{Deserialize, Serialize};

use crate::api::users::model::UserId;

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

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct UpdateTeamSchema {
    pub team_name: Option<String>,
    pub team_description: String
}

#[derive(Debug, Deserialize, Serialize)]
pub struct RegisterDatabaseSchema {
    pub database_name: String,
    pub database_description: String
}

#[derive(Debug, Deserialize, Serialize)]
pub struct RegisterProjectSchema {
    pub project_name: String,
    pub project_description: String
}