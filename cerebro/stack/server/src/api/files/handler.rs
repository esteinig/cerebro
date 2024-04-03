
use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::doc, Collection};
use actix_web::{post, get, web, HttpResponse};

use cerebro_model::api::teams::model::DatabaseId;
use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::schema::RegisterFileSchema;
use cerebro_model::api::users::model::Role;

use crate::api::auth::jwt;
use crate::api::server::AppState;

use crate::api::utils::get_teams_db_collection;
use crate::api::cerebro::handler::get_authorized_database;
use crate::api::files::mongo::get_latest_files_paginated_pipeline;


#[derive(Deserialize)]
struct FileRegisterQuery {  
   // Required for access authorization in user guard middleware
   db: DatabaseId,
}

#[post("/files/register")]
async fn register_file(data: web::Data<AppState>, schema: web::Json<RegisterFileSchema>, query: web::Query<FileRegisterQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to access the file storage", "data": serde_json::json!({})
        }))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, &db.name().to_string(), "files");
    
    match files_collection
    .find_one(doc! { "id": &schema.id }, None)
    .await
    {
        Ok(None) => {},
        Ok(Some(_)) => {
            return HttpResponse::Conflict().json(
                serde_json::json!({"status": "fail", "message": format!("File with the same unique identifier already exists in database"), "data": {}})
            )
        }
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": err.to_string(), "data": {}})
        ),
    }

    let result = files_collection.insert_one(SeaweedFile::from_schema(&schema), None).await;

    match result {
        Ok(_) => HttpResponse::Ok().json(
            serde_json::json!({"status": "success", "message": format!("File registered in team database"), "data": {}})
        ),
        Err(err) => HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string()), "data": {}})
        ),
    }

}


#[derive(Deserialize)]
struct FileListQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    // Paginated files
    page: u32,
    // Limit files returned
    limit: u32,
}


#[get("/files")]
async fn list_files(data: web::Data<AppState>, query: web::Query<FileListQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to access the file storage", "data": serde_json::json!({})
        }))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, &db.name().to_string(), "files");
    
    let pipeline = get_latest_files_paginated_pipeline(query.page as i64, query.limit as i64);
    
    match files_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let files = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match files.is_empty() {
                false => HttpResponse::Ok().json(
                    serde_json::json!({"status": "success", "message": "Files found and listed", "data": files})
                ),
                true => HttpResponse::Ok().json(
                    serde_json::json!({"status": "fail", "message": "Files not found", "data": []})
                )
            }
        },
        Err(err) => HttpResponse::InternalServerError().body(err.to_string()),
    }

}



// Handler configuration
pub fn files_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_file).service(list_files);
}