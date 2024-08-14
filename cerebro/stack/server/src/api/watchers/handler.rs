use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document}, Collection};
use actix_web::{delete, get, post, patch, web, HttpResponse};

use cerebro_model::api::{teams::model::TeamAdminCollection, watchers::{model::ProductionWatcher, response::{DeleteWatcherResponse, ListWatchersResponse, PingWatcherResponse, RegisterWatcherResponse}, schema::RegisterWatcherSchema}};

use crate::api::auth::jwt::{self, JwtDataMiddleware, TeamAccessQuery};
use crate::api::server::AppState;

use crate::api::utils::get_teams_db_collection;
use crate::api::watchers::mongo::get_registered_watchers_pipeline;


#[post("/watcher/register")]
async fn register_watcher(data: web::Data<AppState>, schema: web::Json<RegisterWatcherSchema>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Watchers);
    
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
    // Get a single watcher by identifier for the response
    id: Option<String>
}


#[get("/watcher")]
async fn list_watchers(data: web::Data<AppState>, query: web::Query<WatcherListQuery>, _: web::Query<TeamAccessQuery>, auth_guard: JwtDataMiddleware) -> HttpResponse {
    
    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Watchers);
    
    let watcher = get_registered_watchers_pipeline(&query.id);
    
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
struct WatchersDeleteQuery {  
    // Optional watcher name
    name: Option<String>,
    // Optional watcher location
    location: Option<String>
}


#[delete("/watcher")]
async fn delete_watchers(data: web::Data<AppState>, query: web::Query<WatchersDeleteQuery>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Watchers);

    let mut delete_query = doc! {};

    if let Some(name) = &query.name {
        delete_query.insert("name", name);
    }

    if let Some(location) = &query.location {
        delete_query.insert("location", location);
    }


    match watcher_collection
        .delete_many(delete_query, None) 
        .await
    {
        Ok(delete_result) => {
            if delete_result.deleted_count > 0 {
                HttpResponse::Ok().json(
                    DeleteWatcherResponse::all_deleted()
                )
            } else {
                HttpResponse::NotFound().json(
                    DeleteWatcherResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            DeleteWatcherResponse::server_error(err.to_string())
        )
    }
}

#[delete("/watcher/{id}")]
async fn delete_watcher(data: web::Data<AppState>, id: web::Path<String>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
        
    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Watchers);
    
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


#[patch("/watcher/{id}")]
async fn ping_watcher(data: web::Data<AppState>, id: web::Path<String>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let watcher_collection: Collection<ProductionWatcher> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Watchers);
    
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
    cfg
        .service(register_watcher)
        .service(delete_watchers)
        .service(delete_watcher)
        .service(ping_watcher)
        .service(list_watchers);
}