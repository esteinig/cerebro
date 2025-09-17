use cerebro_model::api::{teams::model::TeamAdminCollection, training::{model::{TrainingPrefetchRecord, TrainingResult, TrainingSessionRecord}, response::{TrainingPrefetchData, TrainingResponse, TrainingSessionData}, schema::{CreateTrainingPrefetch, CreateTrainingSession, PatchTrainingRecord, QueryTrainingData, QueryTrainingPrefetch, TrainingPrefetchOverview}}};

use chrono::{SecondsFormat, Utc};
use futures::stream::TryStreamExt;
use mongodb::{bson::{doc, from_document, to_bson, Bson}, Collection};
use actix_web::{delete, get, patch, post, web, HttpResponse};

use crate::api::{auth::jwt::{self, TeamAccessQuery}, server::AppState, training::gridfs::{delete_from_gridfs, download_prefetch_from_gridfs, find_unique_gridfs_id_by_filename, upload_prefetch_to_gridfs}, utils::get_teams_db_collection};

#[post("/training/prefetch")]
async fn insert_training_data_handler(data: web::Data<AppState>, body: web::Json<CreateTrainingPrefetch>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingData);
    
    // Check for collection and identifier existing already in training database
    if let Ok(Some(_)) = training_collection
        .find_one(doc! {
            "collection": &body.collection,
            "name": &body.name
        })
        .await
    {
        return HttpResponse::Conflict()
            .json(TrainingResponse::<()>::error("A prefetch data entry with this collection and name already exists"));
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
        .json(TrainingResponse::<()>::error(&format!("GridFS upload Failed: {}", err.to_string())))
    };

    let model = TrainingPrefetchRecord::from_request(&body.into_inner());

    match training_collection
        .insert_one(model)
        .await {
            Ok(_) => HttpResponse::Ok().json(TrainingResponse::<()>::created("created training record")),
            Err(err) => {
                // Try to delete the gridfs file if insert fails
                let _ = delete_from_gridfs(&bucket, gridfs_id).await;
                HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Insert Failed: {}", err.to_string())))
            }
        }
}


#[get("/training/prefetch")]
async fn retrieve_training_data_handler(data: web::Data<AppState>, filter: web::Query<QueryTrainingPrefetch>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingData);
    
    let mut query = doc! {};
    if let Some(c) = &filter.collection {
        query.insert("collection", c.clone());
    }
    if let Some(n) = &filter.name {
        query.insert("name", n.clone());
    }
    if let Some(i) = &filter.id {
        query.insert("id", i.clone());
    }

    let mut cursor = match training_collection.find(query).await {
        Ok(cursor) => cursor,
        Err(e) => {
            return HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Query Failed: {}", e)));
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


#[get("/training/prefetch/{session_id}")]
async fn retrieve_training_session_data_handler(data: web::Data<AppState>, session_id: web::Path<String>, filter: web::Query<QueryTrainingData>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_session_collection: Collection<TrainingSessionRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingSessions);
    let training_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingData);
    
    let bucket = data.db.database(&auth_guard.team.admin.database).gridfs_bucket(None);

    // Find the training session
    let session = match training_session_collection
        .find_one(doc! {
            "id": &session_id.into_inner()
        })
        .await
    {
        Ok(Some(session)) => {
            // User guard so that training session are only returned for the correct user
            if session.user_id != auth_guard.user.id {
                return HttpResponse::NotFound().json(TrainingResponse::<()>::not_found("A training session with this identifier does not exist for this user"))
            } else {
                session
            }
        },
        Ok(None) => return HttpResponse::NotFound().json(TrainingResponse::<()>::not_found("A training session with this identifier does not exist")),
        Err(err) => return HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Failed to search for training session: {}", err.to_string())))
    };


    // Convert index and bounds-check
    let idx = match usize::try_from(filter.record) {
        Ok(i) => i,
        Err(_) => {
            return HttpResponse::BadRequest()
                .json(TrainingResponse::<()>::error("Invalid record index"));
        }
    };

    let rec = match session.records.get(idx) {
        Some(r) => r,
        None => {
            return HttpResponse::NotFound()
                .json(TrainingResponse::<()>::not_found("Record index out of range for this session"));
        }
    };
    

    // Fetch the prefetch record by data_id
    let prefetch_record = match training_collection.find_one(doc! { "id": &rec.data_id }).await {
        Ok(Some(prefetch)) => prefetch,
        Ok(None) => {
            return HttpResponse::NotFound().json(TrainingResponse::<()>::not_found(
                "No training data found for the requested record",
            ))
        }
        Err(err) => return HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(
            &format!("Failed to fetch training data: {}", err),
        ))
    };


    let prefetch = match download_prefetch_from_gridfs(&bucket, &prefetch_record.id).await {
        Ok(data) => data,
        Err(err) => {
            return HttpResponse::InternalServerError().json(
                TrainingResponse::<()>::error(&format!("Failed to read GridFS data {}: {}", prefetch_record.id, err.to_string()))
            );
        }
    };
    let data = TrainingPrefetchData::from_query(prefetch_record, prefetch);

    HttpResponse::Ok().json(TrainingResponse::ok((rec, data)))

}


#[get("/training/overview")]
async fn retrieve_training_overview_handler(data: web::Data<AppState>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingData);

    // Pipeline: group, project, sort
    let pipeline = vec![
        doc! {
            "$group": {
                "_id": "$collection",
                "description": { "$first": "$description" },
                "samples": { "$sum": 1 }
            }
        },
        doc! {
            "$project": {
                "_id": 0,
                "collection": "$_id",
                "description": 1,
                "samples": 1
            }
        },
        doc! { "$sort": { "collection": 1 } },
    ];

    let mut cursor = match training_collection.aggregate(pipeline).await {
        Ok(c) => c,
        Err(e) => {
            return HttpResponse::InternalServerError().json(
                TrainingResponse::<()>::error(&format!("Aggregation failed: {}", e)),
            )
        }
    };

    let mut results: Vec<TrainingPrefetchOverview> = Vec::new();
    while let Some(doc) = match cursor.try_next().await {
        Ok(d) => d,
        Err(e) => {
            return HttpResponse::InternalServerError().json(
                TrainingResponse::<()>::error(&format!("Cursor error: {}", e)),
            )
        }
    } {
        match from_document::<TrainingPrefetchOverview>(doc) {
            Ok(item) => results.push(item),
            Err(e) => {
                return HttpResponse::InternalServerError().json(
                    TrainingResponse::<()>::error(&format!("Decode error: {}", e)),
                )
            }
        }
    }

    if results.is_empty() {
        HttpResponse::NotFound()
            .json(TrainingResponse::<()>::not_found("No training prefetch data"))
    } else {
        HttpResponse::Ok().json(TrainingResponse::ok(results))
    }

}


type PrefetchDataId = String;

#[delete("/training/prefetch/{id}")]
async fn delete_training_data_handler(data: web::Data<AppState>, id: web::Path<PrefetchDataId>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingData);

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
            .json(TrainingResponse::<()>::error(&format!("Delete Failed: {}", e)));
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


#[post("/training/session")]
async fn register_training_session_handler(data: web::Data<AppState>, body: web::Json<CreateTrainingSession>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_session_collection: Collection<TrainingSessionRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingSessions);
    let training_data_collection: Collection<TrainingPrefetchRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingData);
    
    // Check for collection data in database and retrieve records (without prefetch data)
    let query = doc! { "collection": &body.collection };
    
    let mut cursor = match training_data_collection.find(query).await {
        Ok(cursor) => cursor,
        Err(e) => {
            return HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Query Failed: {}", e)));
        }
    };

    let bucket = data.db.database(&auth_guard.team.admin.database).gridfs_bucket(None);


    let mut records: Vec<TrainingPrefetchData> = Vec::new();
    while let Some(doc) = cursor.try_next().await.unwrap_or(None) {

        let prefetch = match download_prefetch_from_gridfs(&bucket, &doc.id).await {
            Ok(data) => data,
            Err(err) => {
                return HttpResponse::InternalServerError().json(
                    TrainingResponse::<()>::error(&format!("Failed to read GridFS data {}: {}", doc.id, err.to_string()))
                );
            }
        };

        records.push(TrainingPrefetchData::from_query(doc, prefetch))

    }

    
    if records.is_empty() {
        return HttpResponse::NotFound().json(TrainingResponse::<()>::not_found("No matching training prefetch data for this collection"))
    }

    // Creates the training session record by adding the TrainingRecords for this collection and inserts it into the database
    let model = TrainingSessionRecord::from_request(&body, records, &auth_guard.user.id, &auth_guard.user.name);

    match training_session_collection
        .insert_one(&model)
        .await {
            Ok(_) => HttpResponse::Ok().json(TrainingResponse::<String>::ok(model.id)), // session id is returned
            Err(err) => HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Failed to create training session: {}", err.to_string())))
        }
}


#[get("/training/session/{session_id}")]
async fn retrieve_training_session_handler(data: web::Data<AppState>, session_id: web::Path<String>, _: web::Query<TeamAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let training_records_collection: Collection<TrainingSessionRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingSessions);
    
    // Find the training session
    match training_records_collection
        .find_one(doc! {
            "id": &session_id.into_inner()
        })
        .await
    {
        Ok(Some(session)) => {
            // User guard so that training session are only returned for the correct user
            if session.user_id != auth_guard.user.id {
                HttpResponse::NotFound().json(TrainingResponse::<()>::not_found("A training session with this identifier does not exist for this user"))
            } else {
                HttpResponse::Ok().json(TrainingResponse::<TrainingSessionData>::ok(TrainingSessionData::from_query(session)))
            }
        },
        Ok(None) => HttpResponse::NotFound().json(TrainingResponse::<()>::not_found("A training session with this identifier does not exist")),
        Err(err) => HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Failed to search for training session: {}", err.to_string())))
    }
}

#[patch("/training/session")]
async fn update_training_record_handler(
    data: web::Data<AppState>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
    body: web::Json<PatchTrainingRecord>,
) -> HttpResponse {

    let training_records_collections: Collection<TrainingSessionRecord> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingSessions);

    // Retrieve session and enforce user guard
    let session = match training_records_collections.find_one(doc! { "id": &body.session_id }).await {
        Ok(Some(s)) => s,
        Ok(None) => {
            return HttpResponse::NotFound()
                .json(TrainingResponse::<()>::not_found("A training session with this identifier does not exist"))
        }
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(TrainingResponse::<()>::error(&format!("Failed to search for training session: {}", err)))
        }
    };

    if session.user_id != auth_guard.user.id {
        return HttpResponse::NotFound().json(TrainingResponse::<()>::not_found(
            "A training session with this identifier does not exist for this user",
        ));
    }

    // Build update: always set result; candidates null if None, else set value (including [])
    let mut set_doc = match to_bson(&body.result) {
        Ok(bson_val) => doc! { "records.$.result": bson_val },
        Err(err) => {
            return HttpResponse::InternalServerError().json(
                TrainingResponse::<()>::error(&format!("Failed to encode result: {}", err)),
            )
        }
    };
    
    match &body.candidates {
        Some(cands) => match to_bson(cands) {
            Ok(bson_val) => {
                set_doc.insert("records.$.candidates", bson_val);
            }
            Err(err) => {
                return HttpResponse::InternalServerError().json(
                    TrainingResponse::<()>::error(&format!("Failed to encode candidates: {}", err)),
                )
            }
        },
        None => {
            set_doc.insert("records.$.candidates", Bson::Null);
        }
    }

    set_doc.insert("last_updated", Bson::String(body.record_id.clone()));

    let filter = doc! {
        "id": &body.session_id,
        "records.data_id": &body.record_id,
    };
    let update = doc! { "$set": set_doc };

    match training_records_collections.update_one(filter, update).await {
        Ok(res) if res.matched_count == 0 => HttpResponse::NotFound().json(
            TrainingResponse::<()>::not_found("No matching record id in this session"),
        ),
        Ok(_) => HttpResponse::Ok().json(TrainingResponse::<String>::ok(String::from("Record updated"))),
        Err(err) => HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Failed to update record: {}", err))),
    }
}


#[patch("/training/session/{session_id}")]
async fn complete_training_session_handler(
    data: web::Data<AppState>,
    session_id: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware
) -> HttpResponse {

    let training_records_collections: Collection<TrainingSessionRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingSessions);

    let session_id = session_id.into_inner();

    // Retrieve session and enforce user guard
    let session = match training_records_collections.find_one(doc! { "id": &session_id }).await {
        Ok(Some(s)) => s,
        Ok(None) => {
            return HttpResponse::NotFound()
                .json(TrainingResponse::<()>::not_found("A training session with this identifier does not exist"))
        }
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(TrainingResponse::<()>::error(&format!("Failed to search for training session: {}", err)))
        }
    };

    if session.user_id != auth_guard.user.id {
        return HttpResponse::NotFound().json(TrainingResponse::<()>::not_found(
            "A training session with this identifier does not exist for this user",
        ));
    }

    let filter = doc! {
        "id": &session_id
    };
    let update = doc! { 
        "$set": doc! { "completed": Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true) } 
    };

    match training_records_collections.update_one(filter, update).await {
        Ok(res) if res.matched_count == 0 => return HttpResponse::NotFound().json(
            TrainingResponse::<()>::not_found("No matching record identifier in this session"),
        ),
        Err(err) => return HttpResponse::InternalServerError().json(TrainingResponse::<()>::error(&format!("Failed to update record: {}", err))),
        Ok(_) => {}
    }

    // Compute result and return response
    HttpResponse::Ok().json(TrainingResponse::<TrainingResult>::ok(session.evaluate()))

}


#[get("/training/session/{session_id}/result")]
async fn get_training_session_result_handler(
    data: web::Data<AppState>,
    session_id: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware
) -> HttpResponse {

    let training_records_collections: Collection<TrainingSessionRecord> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::TrainingSessions);

    let session_id = session_id.into_inner();

    // Retrieve session and enforce user guard
    let session = match training_records_collections.find_one(doc! { "id": &session_id }).await {
        Ok(Some(s)) => s,
        Ok(None) => {
            return HttpResponse::NotFound()
                .json(TrainingResponse::<()>::not_found("A training session with this identifier does not exist"))
        }
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(TrainingResponse::<()>::error(&format!("Failed to search for training session: {}", err)))
        }
    };

    if session.user_id != auth_guard.user.id {
        return HttpResponse::NotFound().json(TrainingResponse::<()>::not_found(
            "A training session with this identifier does not exist for this user",
        ));
    }

    // Compute result and return response
    HttpResponse::Ok().json(TrainingResponse::<TrainingResult>::ok(session.evaluate()))

}

pub fn training_config(cfg: &mut web::ServiceConfig) {
    
    cfg.service(insert_training_data_handler)
       .service(retrieve_training_data_handler) 
       .service(delete_training_data_handler)
       .service(retrieve_training_overview_handler)
       .service(register_training_session_handler)
       .service(retrieve_training_session_handler)
       .service(complete_training_session_handler)
       .service(retrieve_training_session_data_handler)
       .service(get_training_session_result_handler)
       .service(update_training_record_handler);

}