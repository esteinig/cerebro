use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document}, Collection};
use actix_web::{delete, get, patch, post, web, HttpResponse};

use cerebro_model::api::{auth::response::UserRoleResponse, pipelines::{model::ProductionPipeline, response::{DeletePipelineResponse, ListPipelinesResponse, PingPipelineResponse, RegisterPipelineResponse}, schema::RegisterPipelineSchema}, teams::model::DatabaseId, utils::HttpMethod};
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
        return HttpResponse::Unauthorized().json(
            UserRoleResponse::unauthorized("/pipeline/register", HttpMethod::Post)
        )
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, 
        Err(error_response) => return error_response
    };

    let pipeline_collection: Collection<ProductionPipeline> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Pipelines);
    
    match pipeline_collection
        .find_one(doc! {
            "$or": [
                { "id": &schema.id },
                { "name": &schema.name, "location": &schema.location }
            ]
        }, None)
        .await
    {
        Ok(None) => {},
        Ok(Some(_)) => return HttpResponse::Conflict().json(
            RegisterPipelineResponse::conflict(&schema.id, &schema.name, &schema.location)
        ),
        Err(err) => return HttpResponse::InternalServerError().json(
            RegisterPipelineResponse::server_error(err.to_string())
        )
    }

    match pipeline_collection
        .insert_one(ProductionPipeline::from_schema(&schema), None)
        .await
    {
        Ok(_) => HttpResponse::Ok().json(
            RegisterPipelineResponse::success(&schema.id)
        ),
        Err(err) => HttpResponse::InternalServerError().json(
            RegisterPipelineResponse::server_error(err.to_string())
        )
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
        return HttpResponse::Unauthorized().json(
            UserRoleResponse::unauthorized("/pipeline", HttpMethod::Get)
        )
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let pipeline_collection: Collection<ProductionPipeline> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Pipelines);
    
    let pipeline = get_registered_pipelines_pipeline();
    
    match pipeline_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            
            let pipelines: Vec<ProductionPipeline> = cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| from_document(doc).ok())
                .collect();

            match pipelines.is_empty() {
                false => HttpResponse::Ok().json(
                    ListPipelinesResponse::success(pipelines)
                ),
                true => HttpResponse::NotFound().json(
                    ListPipelinesResponse::not_found()
                )
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(
            ListPipelinesResponse::server_error(err.to_string())
        )
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
        return HttpResponse::Unauthorized().json(
            UserRoleResponse::unauthorized("/pipeline", HttpMethod::Delete)
        )
    }
    
    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let pipeline_collection: Collection<ProductionPipeline> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Pipelines);
    

    match pipeline_collection
        .find_one_and_delete(
            doc! { "id":  &id.into_inner() },
            None
        )
        .await
    {   
        Ok(deleted) => {
            match deleted {
                Some(pipeline) => HttpResponse::Ok().json(
                    DeletePipelineResponse::success(pipeline)
                ),
                None => HttpResponse::NotFound().json(
                    DeletePipelineResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            DeletePipelineResponse::server_error(err.to_string())
        )
    }
}


#[derive(Deserialize)]
struct PipelinePingQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId,
}

#[patch("/pipeline/{id}")]
async fn ping_pipeline(data: web::Data<AppState>, id: web::Path<String>, query: web::Query<PipelinePingQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(
            UserRoleResponse::unauthorized("/pipeline", HttpMethod::Patch)
        )
    }
    
    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let pipeline_collection: Collection<ProductionPipeline> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Pipelines);

    let new_ping = chrono::Utc::now().to_string();

    match pipeline_collection
        .find_one_and_update(
            doc! { "id":  &id.into_inner()},
            doc! { "$set": { "last_ping": &new_ping } },
            None
        )
        .await
    {   
        Ok(updated) => {
            match updated {
                Some(_) => HttpResponse::Ok().json(
                    PingPipelineResponse::success(new_ping)
                ),
                None => HttpResponse::NotFound().json(
                    PingPipelineResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            PingPipelineResponse::server_error(err.to_string())
        )
    }
}

// Handler configuration
pub fn pipelines_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_pipeline)
       .service(delete_pipeline)
       .service(ping_pipeline)
       .service(list_pipelines);
}