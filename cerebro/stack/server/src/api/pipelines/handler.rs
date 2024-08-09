use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document}, Collection};
use actix_web::{delete, get, post, web, HttpResponse};

use cerebro_model::api::{pipelines::{model::ProductionPipeline, response::{DeletePipelineResponse, ListPipelinesResponse, RegisterPipelineResponse}, schema::RegisterPipelineSchema}, teams::model::DatabaseId};
use cerebro_model::api::users::model::Role;

use crate::api::{auth::jwt, utils::TeamDatabaseInternal};
use crate::api::server::AppState;

use crate::api::utils::get_teams_db_collection;
use crate::api::cerebro::handler::get_authorized_database;
use crate::api::pipelines::mongo::get_registered_pipelines_pipeline;

#[derive(Deserialize)]
struct PipelineRegisterQuery {  
   // Required for access authorization in user guard middleware
   db: DatabaseId,
}

#[post("/pipeline/register")]
async fn register_pipeline(data: web::Data<AppState>, schema: web::Json<RegisterPipelineSchema>, query: web::Query<PipelineRegisterQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to access the file data", "data": serde_json::json!({})
        }))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let files_collection: Collection<ProductionPipeline> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Pipelines);
    
    match files_collection
        .find_one(doc! { "id": &schema.id }, None)
        .await
    {
        Ok(None) => {},
        Ok(Some(_)) => return HttpResponse::Conflict().json(RegisterPipelineResponse::conflict(&schema.id)),
        Err(err) => return HttpResponse::InternalServerError().json(RegisterPipelineResponse::server_error(err.to_string())),
    }

    let result = files_collection.insert_one(ProductionPipeline::from_schema(&schema), None).await;

    match result {
        Ok(_) => HttpResponse::Ok().json(RegisterPipelineResponse::success(&schema.id)),
        Err(err) => HttpResponse::InternalServerError().json(RegisterPipelineResponse::server_error(err.to_string())),
    }

}


#[derive(Deserialize)]
struct PipelineListQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId
}


#[get("/pipeline")]
async fn list_pipelines(data: web::Data<AppState>, query: web::Query<PipelineListQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to access file data", "data": serde_json::json!({})
        }))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let files_collection: Collection<ProductionPipeline> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Pipelines);
    
    let pipeline = get_registered_pipelines_pipeline();
    
    match files_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            
            let files: Vec<ProductionPipeline> = cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| from_document(doc).ok())
                .collect();

            match files.is_empty() {
                false => HttpResponse::Ok().json(ListPipelinesResponse::success(files)),
                true => HttpResponse::NotFound().json(ListPipelinesResponse::not_found())
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(ListPipelinesResponse::server_error(err.to_string()))
    }

}



#[derive(Deserialize)]
struct PipelineDeleteQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId,
}

#[delete("/pipeline/{id}")]
async fn delete_pipeline(data: web::Data<AppState>, id: web::Path<String>, query: web::Query<PipelineDeleteQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to access the file data", "data": serde_json::json!({})
        }))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let files_collection: Collection<ProductionPipeline> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Pipelines);
    
    match files_collection
        .find_one_and_delete(
            doc! { "id":  &id.into_inner()}, None)
        .await
    {   
        Ok(deleted) => {
            match deleted {
                Some(file) => HttpResponse::Ok().json(DeletePipelineResponse::success(file)),
                None => HttpResponse::NotFound().json(DeletePipelineResponse::not_found())
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(DeletePipelineResponse::server_error(err.to_string()))
    }
}



// Handler configuration
pub fn pipelines_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_pipeline)
       .service(delete_pipeline)
       .service(list_pipelines);
}