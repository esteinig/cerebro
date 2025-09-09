use cerebro_model::api::{teams::model::TeamAdminCollection, training::{model::TrainingPrefetchRecord, response::{TrainingPrefetchData, TrainingResponse}, schema::{CreateTrainingPrefetch, QueryTrainingPrefetch}}};

use futures::stream::TryStreamExt;
use mongodb::{bson::doc, Collection};
use actix_web::{get, post, web, HttpResponse, delete};

use crate::api::{auth::jwt::{self, TeamAccessQuery}, server::AppState, training::gridfs::{delete_from_gridfs, download_prefetch_from_gridfs, find_unique_gridfs_id_by_filename, upload_prefetch_to_gridfs}, utils::get_teams_db_collection};

#[post("/training/prefetch")]
async fn insert_training_handler(data: web::Data<AppState>, body: web::Json<CreateTrainingPrefetch>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::Training);
    
    // Check for collection and identifier existing already in training database
    if let Ok(Some(_)) = training_collection
        .find_one(doc! {
            "collection": &body.collection,
            "name": &body.name
        })
        .await
    {
        return HttpResponse::Conflict()
            .json(TrainingResponse::<()>::error("A prefetch with this collection and identifier already exists"));
    }

    let bucket = data.db.database(&auth_guard.team.admin.database).gridfs_bucket(None);

    // Upload the prefetch data to GridFS
    let gridfs_id = match upload_prefetch_to_gridfs(
        &bucket, 
        &body.prefetch, 
        &body.id
    ).await {
        Ok(gridfs_id) => gridfs_id,
        Err(err) => return HttpResponse::InternalServerError()
        .json(TrainingResponse::<()>::error(&format!("GridFS upload failed: {}", err.to_string())))
    };

    let model = TrainingPrefetchRecord::from_request(&body.into_inner());

    match training_collection
        .insert_one(model)
        .await {
            Ok(_) => HttpResponse::Created().json(TrainingResponse::<()>::created("created training record")),
            Err(err) => {
                // Try to delete the gridfs file if insert fails
                let _ = delete_from_gridfs(&bucket, gridfs_id).await;
                HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Insert failed: {}", err.to_string())))
            }
        }
}


#[get("/training/prefetch")]
async fn retrieve_training_handler(data: web::Data<AppState>, filter: web::Query<QueryTrainingPrefetch>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::Training);
    
    let mut query = doc! {};
    if let Some(c) = &filter.collection {
        query.insert("collection", c.clone());
    }
    if let Some(n) = &filter.name {
        query.insert("name", n.clone());
    }

    let mut cursor = match training_collection.find(query).await {
        Ok(cursor) => cursor,
        Err(e) => {
            return HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Query failed: {}", e)));
        }
    };


    let bucket = data.db.database(&auth_guard.team.admin.database).gridfs_bucket(None);
    let mut response_data: Vec<TrainingPrefetchData> = Vec::new();


    while let Some(doc) = cursor.try_next().await.unwrap_or(None) {

        let prefetch = match download_prefetch_from_gridfs(&bucket, &doc.id).await {
            Ok(data) => data,
            Err(err) => {
                return HttpResponse::InternalServerError().json(
                    TrainingResponse::<()>::error(&format!("Failed to read GridFS data {}: {}", doc.id, err.to_string()))
                );
            }
        };

        response_data.push(TrainingPrefetchData::from_query(doc, prefetch))

    }

    if response_data.is_empty() {
        HttpResponse::NotFound().json(TrainingResponse::<()>::not_found("No matching training prefetch data"))
    } else {
        HttpResponse::Ok().json(TrainingResponse::ok(response_data))
    }

}


type PrefetchDataId = String;

#[delete("/training/prefetch/{id}")]
async fn delete_training_handler(data: web::Data<AppState>, id: web::Path<PrefetchDataId>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::Training);

    // Get the document for the GridFS identifier
    let doc = match training_collection
        .find_one(doc! { "id": &id.into_inner() })
        .await
    {   
        Ok(Some(doc)) => doc,
        _ => return HttpResponse::NotFound().json(TrainingResponse::<()>::not_found("A prefetch with this identifier does not exist"))
    };

    // Delete metadata
    if let Err(e) = training_collection.delete_one(doc! { "id": &doc.id }).await {
        return HttpResponse::InternalServerError()
            .json(TrainingResponse::<()>::error(&format!("Delete failed: {}", e)));
    }

    // Delete from GridFS
    let bucket = data.db.database(&auth_guard.team.admin.database).gridfs_bucket(None);

    let gridfs_id = match find_unique_gridfs_id_by_filename(&bucket, &doc.id).await {
        Ok(id) => id,
        Err(err) => return HttpResponse::InternalServerError().json(
            TrainingResponse::<()>::error(&format!("Failed to find unique GridFS identifier {}: {}", doc.id, err.to_string()))
        )
    };

    if let Err(err) = delete_from_gridfs(&bucket, gridfs_id).await {
        return HttpResponse::InternalServerError().json(
            TrainingResponse::<()>::error(&format!("Failed to delete prefetch data from GridFS {}: {}", doc.id, err.to_string()))
        )
    }

    HttpResponse::Ok().json(TrainingResponse::<()>::ok(()))
    
}


pub fn training_config(cfg: &mut web::ServiceConfig) {
    
    cfg.service(insert_training_handler)
       .service(retrieve_training_handler) 
       .service(delete_training_handler);

}