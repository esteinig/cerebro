// Schemas to validate and access request data

use serde::Deserialize;
use crate::api::users::model::Role;

#[derive(Debug, Deserialize)]
pub struct RegisterUserSchema {
    pub name: String,
    pub email: String,
    pub password: String,
    pub verified: bool,
    pub title: Option<String>,
    pub positions: Vec<String>,
    pub roles: Vec<Role>
}

#[derive(Debug, Deserialize)]
pub struct UpdateUserSchema {
    pub name: String,
    pub email: String,
    pub title: Option<String>,
    pub password: Option<String>,
    pub verified: bool,
    pub positions: Vec<String>,
    pub roles: Vec<Role>
}
