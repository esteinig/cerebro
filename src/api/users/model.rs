use chrono::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use crate::api::users::schema::RegisterUserSchema;



#[derive(Error, Debug)]
pub enum UserModelError {
    /// Indicates failure to parse a user authorization role from a str
    #[error("failed to parse user role")]
    ParseUserRole
}

// Singular role authorization access to resources means
// that an `Admin` must also be a `User` - important!
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum Role {
    Bot,
    User,
    Data,
    Admin,
    Report,
}
impl Role {
    pub fn from_str(s: &str) -> Result<Self, UserModelError>  {
        match s {
            "bot" | "Bot" => Ok(Self::Bot),
            "user" | "User" => Ok(Self::User),
            "data" | "Data" => Ok(Self::Data),
            "admin" | "Admin" => Ok(Self::Admin),
            "report" | "Report" => Ok(Self::Report),
            _ => Err(UserModelError::ParseUserRole)
        }
    }
}
impl Default for Role {
    fn default() -> Self {
        Role::User
    }
}


// A type alias to better track user identifier
pub type UserId = String;

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct User {
    pub id: UserId,
    pub name: String,
    pub email: String,
    pub password: String,
    pub title: Option<String>,
    pub roles: Vec<Role>,
    pub positions: Vec<String>,
    pub image: Option<String>,
    pub verified: bool,
    pub created: DateTime<Utc>,
    pub updated: DateTime<Utc>
}
impl User {
    pub fn from(schema: &RegisterUserSchema, password: &str) -> User {
        User {
            id: uuid::Uuid::new_v4().to_string(),
            name: schema.name.to_owned(),
            email: schema.email.to_owned(),
            password: password.to_owned(),
            title: schema.title.to_owned(),
            positions: schema.positions.to_owned(),
            roles: schema.roles.to_owned(),
            image: None,
            verified: schema.verified,
            created: Utc::now(),
            updated: Utc::now()
        }
    }
    // Combined title, name, positions
    pub fn official(&self) -> String {
        let title_str = match &self.title {
            Some(title) => format!("{} ", title),
            None => "".to_string()
        };
        let position_str = match &self.positions.is_empty() {
            true => "".to_string(),
            false => format!(" ({})", self.positions.join(", "))
        };
        format!("{}{}{}", title_str, self.name, position_str)
    }
}

