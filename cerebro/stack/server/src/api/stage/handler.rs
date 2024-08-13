
use cerebro_model::api::{stage::{model::StagedSample, response::{DeleteStagedSampleResponse, ListStagedSamplesResponse, RegisterStagedSampleResponse}}, teams::model::TeamAdminCollection};
use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document}, Collection};
use actix_web::{delete, get, post, web, HttpResponse};

use cerebro_model::api::files::model::SeaweedFile;

use crate::api::auth::jwt::{self, TeamAccessQuery};
use crate::api::server::AppState;

use crate::api::utils::get_teams_db_collection;
use crate::api::stage::mongo::{get_latest_staged_samples_pipeline, get_staged_samples_pipeline};
use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;



#[post("/stage/register")]
async fn register_staged_samples(data: web::Data<AppState>, schema: web::Json<RegisterStagedSampleSchema>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::Files);
    let stage_collection: Collection<StagedSample> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Stage);
    
    let pipeline = get_staged_samples_pipeline(&schema.into_inner());

    let staged_samples: Vec<StagedSample> = match files_collection
    .aggregate(pipeline, None)
    .await
    {
        Ok(cursor) => {
            cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| from_document(doc).ok())
                .collect()
        },
        Err(err) => {
            return HttpResponse::InternalServerError().json(ListStagedSamplesResponse::server_error(err.to_string()));
        }
    };

    let result = stage_collection.insert_many(&staged_samples, None).await;

    match result {
        Ok(_) => {
            let ids = staged_samples.into_iter().map(|x| x.id).collect();
            HttpResponse::Ok().json(RegisterStagedSampleResponse::success(ids))
        },
        Err(err) => HttpResponse::InternalServerError().json(RegisterStagedSampleResponse::server_error(err.to_string())),
    }

}


#[derive(Deserialize)]
struct StageListQuery {  
    // Optional run identifier
    run_id: Option<String>
}


#[get("/stage")]
async fn list_staged_samples(data: web::Data<AppState>, query: web::Query<StageListQuery>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let stage_collection: Collection<StagedSample> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Stage);
    
    let pipeline = get_latest_staged_samples_pipeline(query.run_id.clone());
    
    match stage_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            
            let samples: Vec<StagedSample> = cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| from_document(doc).ok())
                .collect();

            match samples.is_empty() {
                false => HttpResponse::Ok().json(ListStagedSamplesResponse::success(samples)),
                true => HttpResponse::NotFound().json(ListStagedSamplesResponse::not_found())
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(ListStagedSamplesResponse::server_error(err.to_string()))
    }

}


#[delete("/stage/{id}")]
async fn delete_staged_sample(data: web::Data<AppState>, id: web::Path<String>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let files_collection: Collection<StagedSample> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Stage);
    
    match files_collection
        .find_one_and_delete(
            doc! { "id":  &id.into_inner()}, None)
        .await
    {   
        Ok(deleted) => {
            match deleted {
                Some(file) => HttpResponse::Ok().json(DeleteStagedSampleResponse::success(file)),
                None => HttpResponse::NotFound().json(DeleteStagedSampleResponse::not_found())
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(DeleteStagedSampleResponse::server_error(err.to_string()))
    }
}



// Handler configuration
pub fn stage_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_staged_samples)
        .service(delete_staged_sample)
        .service(list_staged_samples);
}