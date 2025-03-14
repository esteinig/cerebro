use actix_web::{web, HttpResponse, HttpRequest};
use cerebro_model::api::teams::model::{Team, TeamAdminCollection};
use cerebro_model::api::utils::AdminCollection;
use mongodb::Collection;

use crate::api::utils::get_teams_db_collection;
use crate::api::{server::AppState, utils::get_cerebro_db_collection};
use cerebro_model::api::logs::model::{RequestLog, LogModule, Action, AccessDetails};

/// A utility function to log to both the admin logs and the team database logs
pub async fn log_database_change(data: &web::Data<AppState>, team: Team, request_log: RequestLog) -> Result<(), HttpResponse> {

    let admin_logs_collection: Collection<RequestLog> = get_cerebro_db_collection(data, AdminCollection::Teams);
    let teams_logs_collection: Collection<RequestLog> = get_teams_db_collection(data, team, TeamAdminCollection::Logs);

    let admin_result = admin_logs_collection.insert_one(request_log.clone()).await;
    match admin_result {
        Ok(_) => {},
        Err(_) => return Err(HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Failed to log into admin database", "data": serde_json::json!({})
        })))
    }

    let team_result = teams_logs_collection.insert_one(request_log).await;
    match team_result {
        Ok(_) => Ok(()),
        Err(_) => Err(HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Failed to log into team database", "data": serde_json::json!({})
        }))),
    }


}


// Login related logs
pub async fn log_user_login(data: &web::Data<AppState>, user_email: &str, action: Action, description: &str, request: &HttpRequest) -> Result<(), mongodb::error::Error> {

    let log = RequestLog::new(
        LogModule::UserLogin,
        action.clone(),
        match action {
            Action::InvalidEmail => true,
            _ => false,
        },
        description.to_string(),
        AccessDetails::new(&request, None, Some(user_email), None, None),
    );

    let logs_collection: Collection<RequestLog> = get_cerebro_db_collection(&data, AdminCollection::Logs);
    match logs_collection.insert_one(log).await {
        Ok(_) => Ok(()),
        Err(err) => Err(err)
    }

}

// Login admin authentication attempts
pub async fn log_admin_auth(data: &web::Data<AppState>, user_id: Option<&str>, email: Option<&str>, action: Action, description: &str, request: &HttpRequest) -> Result<(), mongodb::error::Error> {

    let logs_collection: Collection<RequestLog> = get_cerebro_db_collection(&data, AdminCollection::Logs);

    let log = RequestLog::new(
        LogModule::AdminLogin,
        action,
        true,
        description.to_string(),
        AccessDetails::new(&request, user_id, email, None, None),
    );

    match logs_collection.insert_one(log).await {
        Ok(_) => Ok(()),
        Err(err) => Err(err)
    }

}