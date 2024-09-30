use serde::{Deserialize, Serialize};

use super::model::ProductionTower;

#[derive(Serialize, Deserialize)]
pub struct RegisterTowerResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>
}
impl RegisterTowerResponse {
    pub fn success(id: &str) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Tower registered successfully"),
            data: Some(id.to_string())
        }
    }
    pub fn conflict(id: &str, name: &str, location: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("Tower with identifier '{id}' or name '{name}' and location '{location}' already exists in database"),
            data: Some(id.to_string())
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database operation: {}", error_message),
            data: None
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct ListTowersResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<ProductionTower>>
}
impl ListTowersResponse {
    pub fn success(towers: Vec<ProductionTower>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Tower registrations found in database"),
            data: Some(towers)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No tower registrations found in database"),
            data: None
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct DeleteTowerResponse {
    pub status: String,
    pub message: String,
    pub data: Option<ProductionTower>
}
impl DeleteTowerResponse {
    pub fn success(tower: ProductionTower) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Tower registrations deleted from database"),
            data: Some(tower)
        }
    }
    pub fn all_deleted() -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Tower registrations deleted from database"),
            data: None
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No tower registrations found in database"),
            data: None
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None
        }
    }
}


#[derive(Serialize, Deserialize)]
pub struct PingTowerResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>
}
impl PingTowerResponse {
    pub fn success(last_ping: String) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Tower last ping datetime updated"),
            data: Some(last_ping)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No tower registrations found in database"),
            data: None
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None
        }
    }
}