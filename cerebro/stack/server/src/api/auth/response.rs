use serde::{Deserialize, Serialize};

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

