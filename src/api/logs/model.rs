use actix_web::HttpRequest;
use chrono::{Utc, DateTime};
use serde::{Serialize, Deserialize};

use crate::api::{users::model::UserId, teams::model::{DatabaseId, ProjectId}};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AccessDetails {
    pub user_id: Option<UserId>,
    pub user_email: Option<String>,                 
    pub request_details: RequestDetails,    
    pub database_id: Option<DatabaseId>,  
    pub project_id: Option<ProjectId>,
}
impl AccessDetails {
    pub fn new(request: &HttpRequest, user_id: Option<&str>, user_email: Option<&str>, database_id: Option<&str>, project_id: Option<&str>) -> Self {
        Self {
            request_details: RequestDetails::from_request(&request),
            user_id: user_id.map(String::from),
            user_email: user_email.map(String::from),
            database_id: database_id.map(String::from),
            project_id: project_id.map(String::from),
        }
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RequestDetails {
    pub route: String,
    pub method: String,
    pub headers: Vec<String>
}
impl RequestDetails {
    pub fn from_request(req: &HttpRequest) -> Self {
        Self {
            route: req.path().to_string(),
            method: req.method().to_string(),
            headers: req.headers().iter().map(|header| format!("{}: {}", header.0.to_string(), header.1.to_str().unwrap_or("failed-to-parse").to_string())).collect(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LogModule {
    UserComment,
    UserLogin,
    UserDecision,
    UserAction,
    AdminLogin
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Action {
    // Security actions
    InvalidPassword,
    InvalidEmail,
    InvalidVerification,
    ValidCredentials,
    ValidRole,
    InvalidRole,
    InvalidUser,
    InvalidAdmin,
    InvalidToken,
    // Interface actions
    PriorityTaxonAdded,
    PriorityTaxonRemoved,
    PriorityTaxonModified,
    ReportEntryAdded,
    ReportEntryRemoved,
    // Interface decisions
    PriorityTaxonAccepted,
    PriorityTaxonRejected,
    SampleCommentAdded,
    SampleCommentRemoved,
}


type RequestLogId = String;

/// A generic log entry for logging request actions
/// from authenticated users mainly on data routes
/// Main components are the action and extended 
/// description fields (free to choose) supported
/// by user access and request details
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RequestLog {
    pub id: RequestLogId,
    pub date: DateTime<Utc>,
    pub module: LogModule,
    pub action: Action,
    pub critical: bool,
    pub description: String,
    pub access_details: AccessDetails,
}
impl RequestLog {
    pub fn new(module: LogModule, action: Action, critical: bool, description: String, access_details: AccessDetails) -> Self {
        Self {
            id: uuid::Uuid::new_v4().to_string(),
            date: chrono::Utc::now(),
            module,
            action,
            critical,
            description,
            access_details: access_details.to_owned(),
        }
    }
}

