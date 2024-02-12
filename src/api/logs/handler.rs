use crate::api::{logs::mongo::{get_latest_logs_all_pipeline, get_latest_logs_limit_pipeline}, utils::{get_cerebro_db_collection, get_teams_db_collection}, teams::model::DatabaseId};
use serde::Deserialize;
use mongodb::{bson::doc, Collection};
use crate::api::logs::model::RequestLog;
use crate::api::auth::jwt;
use crate::api::server::AppState;
use actix_web::{get, web, HttpResponse};
use futures::TryStreamExt;
use crate::api::cerebro::handler::get_authorized_database;

#[derive(Deserialize)]
struct AdminLogsQuery {
    limit: Option<i64>,
    critical: Option<bool>
}

#[get("/logs/admin")]
async fn get_admin_logs(data: web::Data<AppState>, query: web::Query<AdminLogsQuery>, _: jwt::JwtAdminMiddleware) -> HttpResponse {

    let pipeline = match query.limit {
        Some(limit) => get_latest_logs_limit_pipeline(limit, match query.critical { Some(v) => v, None => false}),
        None => get_latest_logs_all_pipeline(match query.critical { Some(v) => v, None => false})
    };

    let logs_collection: Collection<RequestLog> = get_cerebro_db_collection(&data, "logs");

    
    match logs_collection
    .aggregate(pipeline, None)
    .await
    {
        Ok(cursor) => {
            let logs = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            HttpResponse::Ok().json(serde_json::json!({"status": "success", "message": "Request logs obtained", "data": { "logs": logs }}))
        },
        Err(err) => HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("Request log query failed: {}", err.to_string())})
        ),
    }
}


#[derive(Deserialize)]
struct TeamLogsQuery {  
    db_id: String,
    limit: Option<i64>,
    critical: Option<bool>
}

#[get("/logs/team")]
async fn get_team_logs(data: web::Data<AppState>, query: web::Query<TeamLogsQuery>, auth_user: jwt::JwtUserMiddleware) -> HttpResponse {

    // Check if user belongs to team with requested name and database identifier to get database name for logs collection

    let db: mongodb::Database = match get_authorized_database(&data, &query.db_id, &auth_user) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let pipeline = match query.limit {
        Some(limit) => get_latest_logs_limit_pipeline(limit, match query.critical { Some(v) => v, None => false}),
        None => get_latest_logs_all_pipeline(match query.critical { Some(v) => v, None => false})
    };

    let logs_collection: Collection<RequestLog> = get_teams_db_collection(&data, &db.name().to_string(), "logs");
    
    match logs_collection
    .aggregate(pipeline, None)
    .await
    {
        Ok(cursor) => {
            let logs = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            HttpResponse::Ok().json(serde_json::json!({"status": "success", "message": "Request logs obtained", "data": { "logs": logs }}))
        },
        Err(err) => HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("Request log query failed: {}", err.to_string())})
        ),
    }
}

// Handler configuration
pub fn logs_config(cfg: &mut web::ServiceConfig) {
    cfg.service(get_admin_logs).service(get_team_logs);
}