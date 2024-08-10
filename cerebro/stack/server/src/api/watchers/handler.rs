use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document}, Collection};
use actix_web::{delete, get, post, patch, web, HttpResponse};

use cerebro_model::api::{auth::response::UserRoleResponse, teams::model::DatabaseId, utils::HttpMethod, watchers::{model::ProductionWatcher, response::{DeleteWatcherResponse, ListWatchersResponse, PingWatcherResponse, RegisterWatcherResponse}, schema::RegisterWatcherSchema}};
use cerebro_model::api::users::model::Role;

use crate::api::{auth::jwt, utils::TeamDatabaseInternal};
use crate::api::server::AppState;

use crate::api::utils::get_teams_db_collection;
use crate::api::cerebro::handler::get_authorized_database;
use crate::api::watchers::mongo::get_registered_watchers_pipeline;

#[derive(Deserialize)]
struct WatcherRegisterQuery {  
   // Required for access authorization in user guard middleware
   db: DatabaseId,
}

#[post("/watcher/register")]
async fn register_watcher(data: web::Data<AppState>, schema: web::Json<RegisterWatcherSchema>, query: web::Query<WatcherRegisterQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(UserRoleResponse::unauthorized("/watcher/register", HttpMethod::Post))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Watchers);
    
    match watcher_collection
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
            RegisterWatcherResponse::conflict(&schema.id, &schema.name, &schema.location)
        ),
        Err(err) => return HttpResponse::InternalServerError().json(
            RegisterWatcherResponse::server_error(err.to_string())
        ),
    }

    match watcher_collection
        .insert_one(ProductionWatcher::from_schema(&schema), None)
        .await 
    {
        Ok(_) => HttpResponse::Ok().json(
            RegisterWatcherResponse::success(&schema.id)
        ),
        Err(err) => HttpResponse::InternalServerError().json(
            RegisterWatcherResponse::server_error(err.to_string())
        ),
    }

}


#[derive(Deserialize)]
struct WatcherListQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId
}


#[get("/watcher")]
async fn list_watchers(data: web::Data<AppState>, query: web::Query<WatcherListQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(UserRoleResponse::unauthorized("/watcher", HttpMethod::Get))
    }

    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Watchers);
    
    let watcher = get_registered_watchers_pipeline();
    
    match watcher_collection
        .aggregate(watcher, None)
        .await
    {
        Ok(cursor) => {
            
            let watchers: Vec<ProductionWatcher> = cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| from_document(doc).ok())
                .collect();

            match watchers.is_empty() {
                false => HttpResponse::Ok().json(
                    ListWatchersResponse::success(watchers)
                ),
                true => HttpResponse::NotFound().json(
                    ListWatchersResponse::not_found()
                )
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(
            ListWatchersResponse::server_error(err.to_string())
        )
    }

}



#[derive(Deserialize)]
struct WatcherDeleteQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId,
}

#[delete("/watcher/{id}")]
async fn delete_watcher(data: web::Data<AppState>, id: web::Path<String>, query: web::Query<WatcherDeleteQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(UserRoleResponse::unauthorized("/watcher", HttpMethod::Delete))
    }
    
    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Watchers);
    
    match watcher_collection
        .find_one_and_delete(
            doc! { "id":  &id.into_inner() }, 
            None
        )
        .await
    {   
        Ok(deleted) => {
            match deleted {
                Some(watcher) => HttpResponse::Ok().json(
                    DeleteWatcherResponse::success(watcher)
                ),
                None => HttpResponse::NotFound().json(
                    DeleteWatcherResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            DeleteWatcherResponse::server_error(err.to_string())
        )
    }
}


#[derive(Deserialize)]
struct WatcherPingQuery {  
    // Required for access authorization in user guard middleware
    db: DatabaseId,
}

#[patch("/watcher/{id}")]
async fn ping_watcher(data: web::Data<AppState>, id: web::Path<String>, query: web::Query<WatcherPingQuery>, auth_guard: jwt::JwtUserMiddleware) -> HttpResponse {
    
    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(UserRoleResponse::unauthorized("/watcher", HttpMethod::Patch))
    }
    
    let db: mongodb::Database = match get_authorized_database(&data, &query.db, &auth_guard) {
        Ok(db) => db, Err(error_response) => return error_response
    };

    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, &db.name().to_string(), TeamDatabaseInternal::Watchers);
    
    let new_ping = chrono::Utc::now().to_string();

    match watcher_collection
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
                    PingWatcherResponse::success(new_ping)
                ),
                None => HttpResponse::NotFound().json(
                    PingWatcherResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            PingWatcherResponse::server_error(err.to_string())
        )
    }
}


// Handler configuration
pub fn watchers_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_watcher)
       .service(delete_watcher)
       .service(ping_watcher)
       .service(list_watchers);
}