//! CIQA server API (Stage 12). Persistence for QC datasets, immutable versioned baselines, and
//! per-run regression results. Mirrors `api::training::handler` one-for-one: team-scoped typed
//! collections via `get_teams_db_collection`, JWT auth guard, `TeamAccessQuery`, GridFS for large
//! payloads with insert rollback, and a uniform `CiqaResponse<T>`. Every mutation writes an
//! append-only audit entry and (for baselines) carries a blake3 integrity hash.

use cerebro_model::api::{
    ciqa::{
        model::{QualityControlBaseline, QualityControlDatasetRecord},
        response::CiqaResponse,
        schema::{CreateCiqaBaseline, CreateCiqaDataset, CreateCiqaRegression, QueryCiqa},
    },
    teams::model::TeamAdminCollection,
};

use actix_web::{delete, get, patch, post, web, HttpResponse};
use chrono::Utc;
use futures::stream::TryStreamExt;
use mongodb::{bson::doc, Collection};

use crate::api::{
    auth::jwt::{self, TeamAccessQuery},
    ciqa::gridfs::{delete_from_gridfs, download_prefetch_from_gridfs, upload_prefetch_to_gridfs},
    server::AppState,
    utils::get_teams_db_collection,
};

/// Best-effort append-only audit of a CIQA mutation (SKILLS.md §4). Failures are logged, not fatal.
async fn write_audit(
    data: &web::Data<AppState>,
    auth: &jwt::JwtDataMiddleware,
    action: &str,
    collection: &str,
    record_id: &str,
) {
    let audit: Collection<mongodb::bson::Document> =
        get_teams_db_collection(data, auth.team.clone(), TeamAdminCollection::AuditLogs);
    let entry = doc! {
        "subsystem": "ciqa",
        "action": action,
        "collection": collection,
        "record_id": record_id,
        "user_id": &auth.user.id,
        "user_name": &auth.user.name,
        "timestamp": Utc::now().to_rfc3339(),
    };
    if let Err(e) = audit.insert_one(entry).await {
        log::warn!("ciqa audit write failed ({action} {collection} {record_id}): {e}");
    }
}

// ── Datasets ─────────────────────────────────────────────────────────────────────────────

#[post("/ciqa/dataset")]
async fn insert_ciqa_dataset_handler(
    data: web::Data<AppState>,
    body: web::Json<CreateCiqaDataset>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let collection: Collection<QualityControlDatasetRecord> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaDatasets);

    if let Ok(Some(_)) = collection
        .find_one(doc! { "dataset": &body.dataset, "version": &body.version, "id": &body.id })
        .await
    {
        return HttpResponse::Conflict().json(CiqaResponse::<()>::conflict(
            "A dataset record with this dataset/version/id already exists",
        ));
    }

    let bucket = data
        .db
        .database(&auth_guard.team.admin.database)
        .gridfs_bucket(None);
    let gridfs_id = match upload_prefetch_to_gridfs(&bucket, &body.prefetch, &body.id).await {
        Ok(id) => id,
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(CiqaResponse::<()>::error(&format!("GridFS upload failed: {}", err)))
        }
    };

    let record = QualityControlDatasetRecord::from_request(&body.into_inner());
    match collection.insert_one(&record).await {
        Ok(_) => {
            write_audit(&data, &auth_guard, "create", "ciqa_datasets", &record.id).await;
            HttpResponse::Ok().json(CiqaResponse::<()>::created("created ciqa dataset record"))
        }
        Err(err) => {
            // roll back the orphaned GridFS blob (training pattern)
            let _ = delete_from_gridfs(&bucket, gridfs_id).await;
            HttpResponse::InternalServerError()
                .json(CiqaResponse::<()>::error(&format!("Insert failed: {}", err)))
        }
    }
}

#[get("/ciqa/dataset")]
async fn list_ciqa_datasets_handler(
    data: web::Data<AppState>,
    filter: web::Query<QueryCiqa>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let collection: Collection<QualityControlDatasetRecord> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaDatasets);

    let mut query = doc! {};
    if let Some(d) = &filter.dataset {
        query.insert("dataset", d.clone());
    }
    if let Some(v) = &filter.version {
        query.insert("version", v.clone());
    }

    let cursor = match collection.find(query).await {
        Ok(c) => c,
        Err(e) => {
            return HttpResponse::InternalServerError()
                .json(CiqaResponse::<()>::error(&format!("Query failed: {}", e)))
        }
    };
    let records: Vec<QualityControlDatasetRecord> = cursor.try_collect().await.unwrap_or_default();
    if records.is_empty() {
        HttpResponse::NotFound().json(CiqaResponse::<()>::not_found("No matching ciqa datasets"))
    } else {
        HttpResponse::Ok().json(CiqaResponse::ok(records))
    }
}

#[get("/ciqa/dataset/{id}")]
async fn get_ciqa_dataset_handler(
    data: web::Data<AppState>,
    path: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let id = path.into_inner();
    let collection: Collection<QualityControlDatasetRecord> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaDatasets);

    let record = match collection.find_one(doc! { "id": &id }).await {
        Ok(Some(r)) => r,
        Ok(None) => {
            return HttpResponse::NotFound()
                .json(CiqaResponse::<()>::not_found("No matching ciqa dataset"))
        }
        Err(e) => {
            return HttpResponse::InternalServerError()
                .json(CiqaResponse::<()>::error(&format!("Query failed: {}", e)))
        }
    };

    let bucket = data
        .db
        .database(&auth_guard.team.admin.database)
        .gridfs_bucket(None);
    match download_prefetch_from_gridfs(&bucket, &record.id).await {
        Ok(prefetch) => HttpResponse::Ok().json(CiqaResponse::ok((record, prefetch))),
        Err(err) => HttpResponse::InternalServerError()
            .json(CiqaResponse::<()>::error(&format!("GridFS read failed: {}", err))),
    }
}

#[delete("/ciqa/dataset/{id}")]
async fn delete_ciqa_dataset_handler(
    data: web::Data<AppState>,
    path: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let id = path.into_inner();
    let collection: Collection<QualityControlDatasetRecord> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaDatasets);

    // Refuse to delete a dataset version that a baseline still references (immutability, D12).
    let baselines: Collection<QualityControlBaseline> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaBaselines);
    if let Ok(Some(record)) = collection.find_one(doc! { "id": &id }).await {
        if let Ok(Some(_)) = baselines
            .find_one(doc! { "dataset": &record.dataset, "dataset_version": &record.version })
            .await
        {
            return HttpResponse::Conflict().json(CiqaResponse::<()>::conflict(
                "Refusing to delete: a baseline still references this dataset version",
            ));
        }
    }

    match collection.delete_one(doc! { "id": &id }).await {
        Ok(res) if res.deleted_count > 0 => {
            let bucket = data
                .db
                .database(&auth_guard.team.admin.database)
                .gridfs_bucket(None);
            if let Ok(file_id) =
                crate::api::ciqa::gridfs::find_unique_gridfs_id_by_filename(&bucket, &id).await
            {
                let _ = delete_from_gridfs(&bucket, file_id).await;
            }
            write_audit(&data, &auth_guard, "delete", "ciqa_datasets", &id).await;
            HttpResponse::Ok().json(CiqaResponse::<()>::created("deleted ciqa dataset record"))
        }
        Ok(_) => HttpResponse::NotFound()
            .json(CiqaResponse::<()>::not_found("No matching ciqa dataset")),
        Err(e) => HttpResponse::InternalServerError()
            .json(CiqaResponse::<()>::error(&format!("Delete failed: {}", e))),
    }
}

// ── Baselines ────────────────────────────────────────────────────────────────────────────

#[post("/ciqa/baseline")]
async fn insert_ciqa_baseline_handler(
    data: web::Data<AppState>,
    body: web::Json<CreateCiqaBaseline>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let collection: Collection<QualityControlBaseline> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaBaselines);

    if let Ok(Some(_)) = collection.find_one(doc! { "id": &body.id }).await {
        return HttpResponse::Conflict()
            .json(CiqaResponse::<()>::conflict("A baseline with this id already exists"));
    }

    // server sets created + content_hash + promoted=false (append-only; D12)
    let baseline = QualityControlBaseline::from_request(&body.into_inner());
    match collection.insert_one(&baseline).await {
        Ok(_) => {
            write_audit(&data, &auth_guard, "create", "ciqa_baselines", &baseline.id).await;
            HttpResponse::Ok().json(CiqaResponse::<()>::created("created ciqa baseline"))
        }
        Err(e) => HttpResponse::InternalServerError()
            .json(CiqaResponse::<()>::error(&format!("Insert failed: {}", e))),
    }
}

#[get("/ciqa/baseline")]
async fn list_ciqa_baselines_handler(
    data: web::Data<AppState>,
    filter: web::Query<QueryCiqa>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let collection: Collection<QualityControlBaseline> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaBaselines);

    let mut query = doc! {};
    if let Some(d) = &filter.dataset {
        query.insert("dataset", d.clone());
    }
    if let Some(v) = &filter.version {
        query.insert("dataset_version", v.clone());
    }

    let cursor = match collection.find(query).await {
        Ok(c) => c,
        Err(e) => {
            return HttpResponse::InternalServerError()
                .json(CiqaResponse::<()>::error(&format!("Query failed: {}", e)))
        }
    };
    let records: Vec<QualityControlBaseline> = cursor.try_collect().await.unwrap_or_default();
    if records.is_empty() {
        HttpResponse::NotFound().json(CiqaResponse::<()>::not_found("No matching ciqa baselines"))
    } else {
        HttpResponse::Ok().json(CiqaResponse::ok(records))
    }
}

// NB: registered before `/ciqa/baseline/{id}` so "active" is not captured as an id.
#[get("/ciqa/baseline/active")]
async fn get_active_ciqa_baseline_handler(
    data: web::Data<AppState>,
    filter: web::Query<QueryCiqa>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let collection: Collection<QualityControlBaseline> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaBaselines);

    let mut query = doc! { "promoted": true };
    if let Some(d) = &filter.dataset {
        query.insert("dataset", d.clone());
    }
    if let Some(v) = &filter.version {
        query.insert("dataset_version", v.clone());
    }

    match collection.find_one(query).await {
        Ok(Some(b)) => HttpResponse::Ok().json(CiqaResponse::ok(b)),
        Ok(None) => HttpResponse::NotFound()
            .json(CiqaResponse::<()>::not_found("No active baseline for this dataset")),
        Err(e) => HttpResponse::InternalServerError()
            .json(CiqaResponse::<()>::error(&format!("Query failed: {}", e))),
    }
}

#[patch("/ciqa/baseline/{id}/promote")]
async fn promote_ciqa_baseline_handler(
    data: web::Data<AppState>,
    path: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let id = path.into_inner();
    let collection: Collection<QualityControlBaseline> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaBaselines);

    let baseline = match collection.find_one(doc! { "id": &id }).await {
        Ok(Some(b)) => b,
        Ok(None) => {
            return HttpResponse::NotFound()
                .json(CiqaResponse::<()>::not_found("No baseline with this id"))
        }
        Err(e) => {
            return HttpResponse::InternalServerError()
                .json(CiqaResponse::<()>::error(&format!("Query failed: {}", e)))
        }
    };

    // Demote any currently-active baseline for this dataset, then promote this one. (For strict
    // atomicity under concurrent promotes, wrap these two updates in a transaction — R2.3.)
    let _ = collection
        .update_many(
            doc! { "dataset": &baseline.dataset, "promoted": true },
            doc! { "$set": { "promoted": false } },
        )
        .await;
    match collection
        .update_one(doc! { "id": &id }, doc! { "$set": { "promoted": true } })
        .await
    {
        Ok(_) => {
            write_audit(&data, &auth_guard, "promote", "ciqa_baselines", &id).await;
            HttpResponse::Ok().json(CiqaResponse::<()>::created("promoted ciqa baseline to active"))
        }
        Err(e) => HttpResponse::InternalServerError()
            .json(CiqaResponse::<()>::error(&format!("Promote failed: {}", e))),
    }
}

#[delete("/ciqa/baseline/{id}")]
async fn delete_ciqa_baseline_handler(
    data: web::Data<AppState>,
    path: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let id = path.into_inner();
    let collection: Collection<QualityControlBaseline> =
        get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::CiqaBaselines);

    match collection.delete_one(doc! { "id": &id }).await {
        Ok(res) if res.deleted_count > 0 => {
            write_audit(&data, &auth_guard, "delete", "ciqa_baselines", &id).await;
            HttpResponse::Ok().json(CiqaResponse::<()>::created("deleted ciqa baseline"))
        }
        Ok(_) => HttpResponse::NotFound()
            .json(CiqaResponse::<()>::not_found("No matching ciqa baseline")),
        Err(e) => HttpResponse::InternalServerError()
            .json(CiqaResponse::<()>::error(&format!("Delete failed: {}", e))),
    }
}

// ── Regression runs ──────────────────────────────────────────────────────────────────────

#[post("/ciqa/regression")]
async fn insert_ciqa_regression_handler(
    data: web::Data<AppState>,
    body: web::Json<CreateCiqaRegression>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let collection: Collection<CreateCiqaRegression> = get_teams_db_collection(
        &data,
        auth_guard.team.clone(),
        TeamAdminCollection::CiqaRegression,
    );

    let mut record = body.into_inner();
    if record.created.is_none() {
        record.created = Some(Utc::now().to_rfc3339());
    }

    match collection.insert_one(&record).await {
        Ok(_) => {
            write_audit(&data, &auth_guard, "create", "ciqa_regression", &record.run_id).await;
            HttpResponse::Ok().json(CiqaResponse::<()>::created("stored ciqa regression result"))
        }
        Err(e) => HttpResponse::InternalServerError()
            .json(CiqaResponse::<()>::error(&format!("Insert failed: {}", e))),
    }
}

#[get("/ciqa/regression/{run_id}")]
async fn get_ciqa_regression_handler(
    data: web::Data<AppState>,
    path: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let run_id = path.into_inner();
    let collection: Collection<CreateCiqaRegression> = get_teams_db_collection(
        &data,
        auth_guard.team.clone(),
        TeamAdminCollection::CiqaRegression,
    );

    match collection.find_one(doc! { "run_id": &run_id }).await {
        Ok(Some(r)) => HttpResponse::Ok().json(CiqaResponse::ok(r)),
        Ok(None) => HttpResponse::NotFound()
            .json(CiqaResponse::<()>::not_found("No regression result for this run_id")),
        Err(e) => HttpResponse::InternalServerError()
            .json(CiqaResponse::<()>::error(&format!("Query failed: {}", e))),
    }
}

/// Register the CIQA handlers. `/ciqa/baseline/active` is registered before `/ciqa/baseline/{id}`
/// so the static segment wins.
pub fn ciqa_config(cfg: &mut web::ServiceConfig) {
    cfg.service(insert_ciqa_dataset_handler)
        .service(list_ciqa_datasets_handler)
        .service(get_ciqa_dataset_handler)
        .service(delete_ciqa_dataset_handler)
        .service(insert_ciqa_baseline_handler)
        .service(list_ciqa_baselines_handler)
        .service(get_active_ciqa_baseline_handler)
        .service(promote_ciqa_baseline_handler)
        .service(delete_ciqa_baseline_handler)
        .service(insert_ciqa_regression_handler)
        .service(get_ciqa_regression_handler);
}
