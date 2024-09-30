
use cerebro_model::api::{towers::model::ProductionTower, stage::{model::StagedSample, response::{DeleteStagedSampleResponse, ListStagedSamplesResponse, RegisterStagedSampleResponse}}, teams::model::TeamAdminCollection};
use serde::Deserialize;
use futures::TryStreamExt;
use mongodb::{bson::{doc, from_document}, Collection};
use actix_web::{delete, get, post, web, HttpResponse};

use cerebro_model::api::files::model::SeaweedFile;

use crate::api::auth::jwt::{self, TeamAccessQuery, TeamProjectAccessQuery};
use crate::api::server::AppState;

use crate::api::utils::{get_teams_db_collection, get_teams_db_stage_collection};
use crate::api::stage::mongo::{get_latest_staged_samples_pipeline, create_staged_samples_pipeline};
use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;

async fn get_tower_from_db(data: &web::Data<AppState>, id: &str, auth_guard: &jwt::JwtDataMiddleware) -> Result<ProductionTower, HttpResponse> {
     // Get the registered pipeline to use the pipeline staging area
     let tower_collection: Collection<ProductionTower> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::Towers);

     match tower_collection
         .find_one(doc! { "id": &id }, None)
         .await
     {
         Ok(Some(tower)) => Ok(tower),
         Ok(None) => Err(HttpResponse::NotFound().json(
             RegisterStagedSampleResponse::tower_not_found(id.to_string())
         )),
         Err(err) => Err(HttpResponse::InternalServerError().json(
             RegisterStagedSampleResponse::server_error(err.to_string())
         ))
     }
}


#[post("/stage/register")]
async fn register_staged_samples(data: web::Data<AppState>, schema: web::Json<RegisterStagedSampleSchema>, access: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let tower = match get_tower_from_db(&data, &schema.tower_id, &auth_guard).await {
        Ok(tower) => tower,
        Err(response) => return response
    };

    // Get the registered files from the staged sample request and transform into StagedSample
    let aggregate_pipeline = create_staged_samples_pipeline(&schema.into_inner(), &tower, &access.db, &access.project);
    let files_collection: Collection<SeaweedFile> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::Files);

    let staged_samples: Vec<StagedSample> = match files_collection
        .aggregate(aggregate_pipeline, None)
        .await
    {
        Ok(cursor) => {
            cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| {
                    let staged_sample: Option<StagedSample> = from_document(doc).ok();
                    if let Some(mut staged_sample) = staged_sample {
                        staged_sample.id = uuid::Uuid::new_v4().to_string(); // assign a unique staging identifier!
                        Some(staged_sample)
                    } else {
                        None
                    }
                    
                })
                .collect::<Vec<StagedSample>>()
        },
        Err(err) => {
            return HttpResponse::InternalServerError().json(
                RegisterStagedSampleResponse::server_error(err.to_string())
            );
        }
    };

    log::info!("{:#?}", staged_samples);

    // Insert the StagedSample into the tower stage collection
    let stage_collection: Collection<StagedSample> = get_teams_db_stage_collection(&data, auth_guard.team, &tower.stage);

    match stage_collection
        .insert_many(&staged_samples, None)
        .await
    {
        Ok(_) => {
            let ids = staged_samples.into_iter().map(|x| x.id).collect();
            HttpResponse::Ok().json(
                RegisterStagedSampleResponse::success(ids)
            )
        },
        Err(err) => HttpResponse::InternalServerError().json(
            RegisterStagedSampleResponse::server_error(err.to_string())
        ),
    }

}


#[derive(Deserialize)]
struct StageListQuery {  
    // Optional run identifier
    run_id: Option<String>,
    // Optional sample identifier
    sample_id: Option<String>
}


#[get("/stage/{tower_id}")]
async fn list_staged_samples(data: web::Data<AppState>, tower_id: web::Path<String>, query: web::Query<StageListQuery>,  _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let tower = match get_tower_from_db(&data, &tower_id, &auth_guard).await {
        Ok(tower) => tower,
        Err(response) => return response
    };

    let stage_collection: Collection<StagedSample> = get_teams_db_stage_collection(&data, auth_guard.team, &&tower.stage);
    
    let aggregate_pipeline = get_latest_staged_samples_pipeline(query.run_id.clone(), query.sample_id.clone());
    
    match stage_collection
        .aggregate(aggregate_pipeline, None)
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
                false => HttpResponse::Ok().json(
                    ListStagedSamplesResponse::success(samples)
                ),
                true => HttpResponse::NotFound().json(
                    ListStagedSamplesResponse::not_found()
                )
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(ListStagedSamplesResponse::server_error(err.to_string()))
    }
}


#[derive(Deserialize)]
struct StageDeleteQuery {  
    // Optional staged sample identifier
    stage_id: Option<String>,
    // Optional run identifier
    run_id: Option<String>,
    // Optional sample identifier
    sample_id: Option<String>
}

#[delete("/stage/{tower_id}")]
async fn delete_staged_samples(data: web::Data<AppState>, tower_id: web::Path<String>, query: web::Query<StageDeleteQuery>,   _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    

    let tower = match get_tower_from_db(&data, &tower_id, &auth_guard).await {
        Ok(tower) => tower,
        Err(response) => return response
    };


    let stage_collection: Collection<StagedSample> = get_teams_db_stage_collection(&data, auth_guard.team, &tower.stage);
    
    let mut delete_query = doc! {};


    if let Some(id) = &query.stage_id {
        delete_query.insert("id", id);
    }

    if let Some(name) = &query.run_id {
        delete_query.insert("run_id", name);
    }
    
    if let Some(location) = &query.sample_id {
        delete_query.insert("sample_id", location);
    }


    match stage_collection
        .delete_many(delete_query, None) 
        .await
    {
        Ok(delete_result) => {
            if delete_result.deleted_count > 0 {
                HttpResponse::Ok().json(
                    DeleteStagedSampleResponse::all_deleted()
                )
            } else {
                HttpResponse::NotFound().json(
                    DeleteStagedSampleResponse::not_found()
                )
            }
        }
        Err(err) => HttpResponse::InternalServerError().json(
            DeleteStagedSampleResponse::server_error(err.to_string())
        )
    }
}


// Handler configuration
pub fn stage_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_staged_samples)
        .service(delete_staged_samples)
        .service(list_staged_samples);
}