use serde::{Deserialize, Serialize};

use super::model::ProductionWatcher;

#[derive(Serialize, Deserialize)]
pub struct RegisterWatcherResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>
}
impl RegisterWatcherResponse {
    pub fn success(id: &str) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Watcher registered successfully"),
            data: Some(id.to_string())
        }
    }
    pub fn conflict(id: &str, name: &str, location: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("Watcher with identifier '{id}' or name '{name}' and location '{location}' already exists in database"),
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
pub struct ListWatchersResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<ProductionWatcher>>
}
impl ListWatchersResponse {
    pub fn success(watchers: Vec<ProductionWatcher>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Watcher registrations found in database"),
            data: Some(watchers)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No watcher registrations found in database"),
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
pub struct DeleteWatcherResponse {
    pub status: String,
    pub message: String,
    pub data: Option<ProductionWatcher>
}
impl DeleteWatcherResponse {
    pub fn success(watchers: ProductionWatcher) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Watcher registrations deleted from database"),
            data: Some(watchers)
        }
    }
    pub fn all_deleted() -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Watcher registrations deleted from database"),
            data: None
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No watcher registrations found in database"),
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
pub struct PingWatcherResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>
}
impl PingWatcherResponse {
    pub fn success(last_ping: String) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Watcher registrations deleted from database"),
            data: Some(last_ping)
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No watcher registrations found in database"),
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