use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct TokenDetails {
    pub token: Option<String>,
    pub token_uuid: uuid::Uuid,
    pub user_id: uuid::Uuid,
    pub expires_in: Option<i64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct TokenClaims {
    pub token_uuid: String,
    pub sub: String,
    pub exp: i64,
    pub iat: i64,
    pub nbf: i64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct AuthLoginSchema {
    pub email: String,
    pub password: String,
}

#[derive(Debug, Deserialize)]
pub struct AuthOneTimeCheckSchema {
    pub access_token: String,
}


#[derive(Debug, Deserialize)]
pub struct AuthResetPasswordSchema {
    pub access_token: String,
    pub password: String,
}
