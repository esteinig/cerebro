
use cerebro_model::api::{files::{response::{DeleteFileResponse, DeleteFilesResponse, ListFilesResponse, RegisterFileResponse, UpdateTagsResponse}, schema::UpdateFileTagsSchema}, teams::model::TeamAdminCollection};
use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document, to_bson}, Collection};
use actix_web::{delete, get, patch, post, web, HttpResponse};

use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::schema::RegisterFileSchema;

use crate::api::auth::jwt::{self, TeamAccessQuery};
use crate::api::server::AppState;

use crate::api::utils::get_teams_db_collection;
use crate::api::files::mongo::get_latest_files_paginated_pipeline;


#[post("/files/register")]
async fn register_file(data: web::Data<AppState>, schema: web::Json<RegisterFileSchema>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);
    
    match files_collection
        .find_one(doc! { "id": &schema.id })
        .await
    {
        Ok(None) => {},
        Ok(Some(_)) => return HttpResponse::Conflict().json(RegisterFileResponse::conflict(&schema.id)),
        Err(err) => return HttpResponse::InternalServerError().json(RegisterFileResponse::server_error(err.to_string())),
    }

    let result = files_collection.insert_one(SeaweedFile::from_schema(&schema)).await;

    match result {
        Ok(_) => HttpResponse::Ok().json(RegisterFileResponse::success(&schema.id)),
        Err(err) => HttpResponse::InternalServerError().json(RegisterFileResponse::server_error(err.to_string())),
    }
}


#[derive(Deserialize)]
struct FileListQuery {  
    // Paginated files
    page: u32,
    // Limit files returned - can be zero for all files
    limit: u32,
    // Optional run identifier
    run_id: Option<String>,
    // Optional watcher identifier
    watcher_id: Option<String>,
   
}


#[get("/files")]
async fn list_files(data: web::Data<AppState>, query: web::Query<FileListQuery>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);
    
    let pipeline = get_latest_files_paginated_pipeline(query.run_id.clone(), query.watcher_id.clone(), query.page as i64, query.limit as i64);
    
    match files_collection
        .aggregate(pipeline)
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


#[patch("/files/tags")]
async fn update_tags(data: web::Data<AppState>, schema: web::Json<UpdateFileTagsSchema>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);
    
    match files_collection
        .update_many(
            doc! { "id": { "$in": &schema.ids } },
            doc! { "$set": { "tags": to_bson(&schema.tags).unwrap() } }
        )
        .await
    {   
        Ok(update_result) => {
            if update_result.matched_count > 0 {
                HttpResponse::Ok().json(UpdateTagsResponse::success())
            } else {
                HttpResponse::NotFound().json(UpdateTagsResponse::not_found())
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            UpdateTagsResponse::server_error(err.to_string())
        )
    }
}



#[derive(Deserialize)]
struct FilesDeleteQuery {  
    // Optional run identifier
    run_id: Option<String>,
    // Optional sample identifier
    sample_id: Option<String>,
    // Delete all option
    all: Option<bool>
}

#[delete("/files")]
async fn delete_files(
    data: web::Data<AppState>, query: web::Query<FilesDeleteQuery>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {

    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    let mut delete_query = doc! {}; // Default to deleting all if `all` is true

    if let Some(true) = query.all {
        // If `all` is true, delete all documents
        delete_query = doc! {};
    } else {

        if query.run_id.is_none() && query.sample_id.is_none() {
            return HttpResponse::BadRequest().json(DeleteFilesResponse::invalid_query());
        }

        // Build query based on provided parameters
        if let Some(run_id) = &query.run_id {
            delete_query.insert("run_id", run_id);
        }

        if let Some(sample_id) = &query.sample_id {
            delete_query.insert("sample_id", sample_id);
        }
    }

    // Retrieve IDs of files to delete
    match files_collection
        .find(delete_query.clone())
        .await
    {
        Ok(cursor) => {
            let deleted_ids = cursor.try_collect().await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .map(|f| f.fid)
                .collect();

            // Proceed with deletion
            match files_collection.delete_many(delete_query).await {
                Ok(_) => {
                    HttpResponse::Ok().json(DeleteFilesResponse::success(deleted_ids))
                }
                Err(err) => HttpResponse::InternalServerError().json(
                    DeleteFilesResponse::server_error(err.to_string()),
                ),
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            DeleteFilesResponse::server_error(err.to_string()),
        ),
    }
}

#[delete("/files/{id}")]
async fn delete_file(data: web::Data<AppState>, id: web::Path<String>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);
    
    match files_collection
        .find_one_and_delete(
            doc! { "id":  &id.into_inner()}
        )
        .await
    {   
        Ok(deleted) => {
            match deleted {
                Some(file) => HttpResponse::Ok().json(
                    DeleteFileResponse::success(file)
                ),
                None => HttpResponse::NotFound().json(
                    DeleteFileResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            DeleteFileResponse::server_error(err.to_string())
        )
    }
}

// Handler configuration
pub fn files_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_file)
        .service(delete_file)
        .service(delete_files)
        .service(update_tags)
        .service(list_files);
}