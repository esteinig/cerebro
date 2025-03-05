use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document}, Collection};
use actix_web::{delete, get, patch, post, web, HttpResponse};

use cerebro_model::api::{towers::{model::ProductionTower, response::{DeleteTowerResponse, ListTowersResponse, PingTowerResponse, RegisterTowerResponse}, schema::RegisterTowerSchema}, teams::model::TeamAdminCollection};
use crate::api::auth::jwt::{self, TeamAccessQuery};
use crate::api::server::AppState;

use crate::api::utils::get_teams_db_collection;
use crate::api::towers::mongo::get_registered_towers_pipeline;


#[post("/tower/register")]
async fn register_tower(data: web::Data<AppState>, schema: web::Json<RegisterTowerSchema>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let tower_collection: Collection<ProductionTower> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Towers);
    
    match tower_collection
        .find_one(doc! {
            "$or": [
                { "id": &schema.id },
                { "name": &schema.name, "location": &schema.location }
            ]
        })
        .await
    {
        Ok(None) => {},
        Ok(Some(_)) => return HttpResponse::Conflict().json(
            RegisterTowerResponse::conflict(&schema.id, &schema.name, &schema.location)
        ),
        Err(err) => return HttpResponse::InternalServerError().json(
            RegisterTowerResponse::server_error(err.to_string())
        )
    }

    match tower_collection
        .insert_one(ProductionTower::from_schema(&schema))
        .await
    {
        Ok(_) => HttpResponse::Ok().json(
            RegisterTowerResponse::success(&schema.id)
        ),
        Err(err) => HttpResponse::InternalServerError().json(
            RegisterTowerResponse::server_error(err.to_string())
        )
    }

}


#[derive(Deserialize)]
struct TowerListQuery {  
    // Get a single tower by identifier for the response
    id: Option<String>
}


#[get("/tower")]
async fn list_tower(data: web::Data<AppState>, query: web::Query<TowerListQuery>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let tower_collection: Collection<ProductionTower> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Towers);
    
    let tower = get_registered_towers_pipeline(&query.id);
    
    match tower_collection
        .aggregate(tower)
        .await
    {
        Ok(cursor) => {
            
            let towers: Vec<ProductionTower> = cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| from_document(doc).ok())
                .collect();

            match towers.is_empty() {
                false => HttpResponse::Ok().json(
                    ListTowersResponse::success(towers)
                ),
                true => HttpResponse::NotFound().json(
                    ListTowersResponse::not_found()
                )
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(
            ListTowersResponse::server_error(err.to_string())
        )
    }

}


#[derive(Deserialize)]
struct TowerDeleteQuery {  
    // Optional tower name
    name: Option<String>,
    // Optional tower location
    location: Option<String>
}


#[delete("/tower")]
async fn delete_towers(data: web::Data<AppState>, query: web::Query<TowerDeleteQuery>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let tower_collection: Collection<ProductionTower> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Towers);
    
    let mut delete_query = doc! {};

    if let Some(name) = &query.name {
        delete_query.insert("name", name);
    }
    
    if let Some(location) = &query.location {
        delete_query.insert("location", location);
    }

    match tower_collection
        .delete_many(delete_query) 
        .await
    {
        Ok(delete_result) => {
            if delete_result.deleted_count > 0 {
                HttpResponse::Ok().json(DeleteTowerResponse::all_deleted())
            } else {
                HttpResponse::NotFound().json(DeleteTowerResponse::not_found())
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(DeleteTowerResponse::server_error(err.to_string()))
    }
}


#[delete("/tower/{id}")]
async fn delete_tower(data: web::Data<AppState>, id: web::Path<String>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let tower_collection: Collection<ProductionTower> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Towers);

    match tower_collection
        .find_one_and_delete(
            doc! { "id":  &id.into_inner() }
        )
        .await
    {   
        Ok(deleted) => {
            match deleted {
                Some(tower) => HttpResponse::Ok().json(
                    DeleteTowerResponse::success(tower)
                ),
                None => HttpResponse::NotFound().json(
                    DeleteTowerResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            DeleteTowerResponse::server_error(err.to_string())
        )
    }
}


#[patch("/tower/{id}")]
async fn ping_tower(data: web::Data<AppState>, id: web::Path<String>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let tower_collection: Collection<ProductionTower> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Towers);

    let new_ping = chrono::Utc::now().to_string();

    match tower_collection
        .find_one_and_update(
            doc! { "id":  &id.into_inner()},
            doc! { "$set": { "last_ping": &new_ping } }
        )
        .await
    {   
        Ok(updated) => {
            match updated {
                Some(_) => HttpResponse::Ok().json(
                    PingTowerResponse::success(new_ping)
                ),
                None => HttpResponse::NotFound().json(
                    PingTowerResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            PingTowerResponse::server_error(err.to_string())
        )
    }
}

// Handler configuration
pub fn towers_config(cfg: &mut web::ServiceConfig) {
    cfg
        .service(register_tower)
        .service(delete_tower)
        .service(delete_towers)
        .service(ping_tower)
        .service(list_tower);
}