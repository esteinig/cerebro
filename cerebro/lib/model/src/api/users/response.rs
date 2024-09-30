use chrono::prelude::*;
use serde::{Serialize, Deserialize};
use crate::api::{users::model::{UserId, Role}, teams::model::Team};

#[derive(Debug, Deserialize, Serialize)]
pub struct FilteredUser {
    pub id: UserId,
    pub name: String,
    pub email: String,
    pub title: Option<String>,
    pub roles: Vec<Role>,
    pub positions: Vec<String>,
    pub image: Option<String>,
    pub verified: bool,
    pub created: DateTime<Utc>,
    pub updated: DateTime<Utc>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct UserData {
    pub user: FilteredUser,
}


#[derive(Serialize, Deserialize, Debug)]
pub struct TeamsData {
    pub teams: Vec<Team>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TeamData {
    pub team: Team,
}


#[derive(Serialize, Deserialize, Debug)]
pub struct UserSelfResponse {
    pub status: String,
    pub message: String,
    pub data: UserData,
}


#[derive(Serialize, Deserialize, Debug)]
pub struct UserSelfTeamsResponse {
    pub status: String,
    pub message: String,
    pub data: TeamsData,
}


#[derive(Serialize, Deserialize, Debug)]
pub struct UserSelfTeamResponse {
    pub status: String,
    pub message: String,
    pub data: TeamData,
}