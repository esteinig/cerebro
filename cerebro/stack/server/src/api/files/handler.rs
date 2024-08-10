
use cerebro_model::api::files::response::{DeleteFileResponse, ListFilesResponse, RegisterFileResponse};
use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document}, Collection};
use actix_web::{delete, get, post, web, HttpResponse};

use cerebro_model::api::teams::model::DatabaseId;
use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::schema::RegisterFileSchema;
use cerebro_model::api::users::model::Role;

use crate::api::{auth::jwt, utils::TeamDatabaseInternal};
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
            "status": "fail", "message": "You do not have permission to access the file data", "data": serde_json::json!({})
        }))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Files);
    
    match files_collection
        .find_one(doc! { "id": &schema.id }, None)
        .await
    {
        Ok(None) => {},
        Ok(Some(_)) => return HttpResponse::Conflict().json(RegisterFileResponse::conflict(&schema.id)),
        Err(err) => return HttpResponse::InternalServerError().json(RegisterFileResponse::server_error(err.to_string())),
    }

    let result = files_collection.insert_one(SeaweedFile::from_schema(&schema), None).await;

    match result {
        Ok(_) => HttpResponse::Ok().json(RegisterFileResponse::success(&schema.id)),
        Err(err) => HttpResponse::InternalServerError().json(RegisterFileResponse::server_error(err.to_string())),
    }

}


#[derive(Deserialize)]
struct FileListQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    // Optional run identifier
    run_id: Option<String>,
    // Optional watcher identifier
    watcher_id: Option<String>,
    // Paginated files
    page: u32,
    // Limit files returned
    limit: u32,
}


#[get("/files")]
async fn list_files(data: web::Data<AppState>, query: web::Query<FileListQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to access file data", "data": serde_json::json!({})
        }))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Files);
    
    let pipeline = get_latest_files_paginated_pipeline(query.run_id.clone(), query.watcher_id.clone(), query.page as i64, query.limit as i64);
    
    match files_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            
            let files: Vec<SeaweedFile> = cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| from_document(doc).ok())
                .collect();

            match files.is_empty() {
                false => HttpResponse::Ok().json(ListFilesResponse::success(files)),
                true => HttpResponse::NotFound().json(ListFilesResponse::not_found())
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(ListFilesResponse::server_error(err.to_string()))
    }

}



#[derive(Deserialize)]
struct FilesDeleteQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId,
}

#[delete("/files/{id}")]
async fn delete_file(data: web::Data<AppState>, id: web::Path<String>, query: web::Query<FilesDeleteQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to access the file data", "data": serde_json::json!({})
        }))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Files);
    
    match files_collection
        .find_one_and_delete(
            doc! { "id":  &id.into_inner()}, None)
        .await
    {   
        Ok(deleted) => {
            match deleted {
                Some(file) => HttpResponse::Ok().json(DeleteFileResponse::success(file)),
                None => HttpResponse::NotFound().json(DeleteFileResponse::not_found())
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(DeleteFileResponse::server_error(err.to_string()))
    }
}



// Handler configuration
pub fn files_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_file)
       .service(delete_file)
       .service(list_files);
}