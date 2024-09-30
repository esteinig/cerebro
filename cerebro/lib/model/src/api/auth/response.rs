use serde::{Deserialize, Serialize};

use crate::api::utils::HttpMethod;

#[derive(Debug, Serialize, Deserialize)]
pub struct AuthLoginResponseSuccess {
    pub access_token: String,
    pub refresh_token: String,
    pub status: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct AuthRefreshResponseSuccess {
    pub access_token: String,
    pub status: String,
}


#[derive(Serialize, Deserialize)]
pub struct UserRoleResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>
}
impl UserRoleResponse {
    pub fn unauthorized(route: &str, method: HttpMethod) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("You do not have permission to access {route} ({method})"),
            data: None
        }
    }
}