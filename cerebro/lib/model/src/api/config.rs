use std::path::PathBuf;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum ConfigError {
    #[error("failed to open config file")]
    ConfigFileInputInvalid,
    #[error("failed to create config file")]
    ConfigFileOutputInvalid,
    #[error("failed to deserialize config file ({0})")]
    ConfigFileNotDeserialized(#[source] toml::de::Error),
    #[error("failed to serialize config file ({0})")]
    ConfigFileNotSerialized(#[source] toml::ser::Error),
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseConfig {
    pub connections: DatabaseConnectionConfig,
    pub names: DatabaseNameConfig
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseConnectionConfig {
    pub mongodb_uri: String,
    pub redis_auth_session_uri: String,
    pub redis_auth_onetime_uri: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseNameConfig {
    // Administrative database collections 
    pub admin_database_name: String,
    pub admin_database_user_collection: String,
    pub admin_database_team_collection: String,
    pub admin_database_logs_collection: String,
    pub admin_database_scheduler_collection: String,
    pub admin_database_jobs_collection: String,
    pub admin_database_locks_collection: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TokenConfig {
    pub encryption: TokenEncryptionConfig,
    pub expiration: TokenExpirationConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TokenExpirationConfig {
    pub access_max_age: i64,
    pub access_max_age_bot: i64,
    pub refresh_max_age: i64,
    pub refresh_max_age_bot: i64,
    pub onetime_max_age_email: i64,
    pub onetime_max_age_password: i64,
    pub onetime_max_age_reset: i64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TokenEncryptionConfig {
    pub access_private_key: String,
    pub access_public_key: String,
    pub refresh_private_key: String,
    pub refresh_public_key: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SecurityConfig { 
    pub components: SecurityComponentsConfig,
    pub token: TokenConfig,
    pub cors: SecurityCorsConfig,
    pub cookies: SecurityCookiesConfig
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SecurityComponentsConfig { 
    pub email: bool,            // verification and password emails
    pub comments: bool,         // sample comments, priority taxon submission and decision comments
    pub annotation: bool,       // sample annotations including fields in sample sheet parser
    pub report_header: bool,    // sample or patient data in report header
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SecurityCookiesConfig { 
    pub domain: String,                  // Cookie domain settings for sharing between subdomains
    pub secure: bool                     // Cookie secure attribute set; requests through SSL/TLS 
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SecurityCorsConfig {
    // Cross-origin access from application; the public
    // application origin is also inserted into emails!
    pub app_origin_public_url: String,
    pub app_origin_docker_url: String,
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SmtpConfig {
    // Email server
    pub host: String,
    pub port: u16,
    pub username: String,
    pub password: String,
    pub from: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub database: DatabaseConfig,
    pub security: SecurityConfig,
    pub smtp: SmtpConfig,
}

impl Config {
    pub fn from_toml(path: &PathBuf) -> Result<Self, ConfigError> {
        let toml_str = std::fs::read_to_string(path).map_err(|_| ConfigError::ConfigFileInputInvalid)?;
        let config = toml::from_str(&toml_str).map_err(|err| ConfigError::ConfigFileNotDeserialized(err))?;
        Ok(config)
    }
    pub fn to_toml(&self, path: &PathBuf) -> Result<(), ConfigError> {
        let config = toml::to_string_pretty(&self).map_err(|err| ConfigError::ConfigFileNotSerialized(err))?;
        std::fs::write(path, config).map_err(|_| ConfigError::ConfigFileOutputInvalid)
    }
}