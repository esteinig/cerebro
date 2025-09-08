use std::time::Duration;

use actix_web::{post, get, web, HttpResponse};
use cerebro_model::api::jobs::schema::{EnqueueJobRequest, ScheduleJobRequest};
use mongodb::bson::{doc, Bson, Uuid};
use mongodb::Collection;
use futures::TryStreamExt;
use crate::api::auth::jwt::JwtAdminMiddleware;
use crate::api::server::AppState;
use crate::api::utils::get_cerebro_db_collection;
use cerebro_model::api::utils::AdminCollection;

use cerebro_model::api::jobs::model::{JobRunDoc, ScheduleJob};
use cerebro_model::api::jobs::response::{
    CreateScheduleResponse, EnqueueJobResponse, JobCompletionResponse, JobsStatusData, JobsStatusResponse, ScheduleJobSummary
};

/// Build a Faktory job whose args is a one-element JSON array
fn job_with_val(kind: &str, val: serde_json::Value) -> faktory::Job {
    faktory::Job::new(kind, vec![val])
}


#[post("/jobs/enqueue")]
pub async fn enqueue_now_handler(
    data: web::Data<AppState>,
    body: web::Json<EnqueueJobRequest>,
    _: JwtAdminMiddleware,
) -> HttpResponse {
    
    let b = body.into_inner();

    // Build Faktory client
    let mut c = match faktory::Client::connect().await {
        Ok(c) => c,
        Err(e) => {
            return HttpResponse::InternalServerError()
                .json(EnqueueJobResponse::server_error(e.to_string()))
        }
    };

    let queue = b.queue.unwrap_or("default".to_string());

    // Build the job and capture JID *before* enqueue
    let mut job = job_with_val(&b.kind, b.args).on_queue(&queue);
    let jid = job.id().to_string();

    job.retry = b.retry.map(|x| x as isize);
    job.reserve_for = b.reserve_for_seconds.map(|x| Duration::from_secs(x as u64));

    // Write an initial JobRunDoc (status = queued)
    let runs: Collection<JobRunDoc> = get_cerebro_db_collection(&data, AdminCollection::Jobs);
    let now = chrono::Utc::now();
    let id = uuid::Uuid::new_v4();

    let doc = JobRunDoc {
        id: Some(id),                 
        jid: jid.clone(),
        kind: b.kind.clone(),
        queue: queue.clone(),
        status: "queued".into(),
        result: None,
        error: None,
        created_at: now,
        updated_at: now,
    };

    if let Err(e) = runs.insert_one(&doc).await {
        return HttpResponse::InternalServerError()
            .json(EnqueueJobResponse::server_error(format!("jobrun insert: {e}")));
    }

    // Enqueue to Faktory
    if let Err(e) = c.enqueue(job).await {
        // Best-effort: mark failed-to-enqueue
        let _ = runs
            .update_one(
                doc! { "jid": &jid },
                doc! { "$set": { "status": "failed", "error": format!("enqueue: {e}"), "updated_at": chrono::Utc::now() } }
            )
            .await;

        return HttpResponse::InternalServerError()
            .json(EnqueueJobResponse::server_error(e.to_string()));
    }

    // Return id [uuid string]
    HttpResponse::Accepted().json(EnqueueJobResponse::success(&id.to_string()))
}

#[get("/jobs/enqueue/{id}")]
pub async fn job_completion_handler(
    data: web::Data<AppState>,
    path: web::Path<String>,
    _: JwtAdminMiddleware,
) -> HttpResponse {
    let id = path.into_inner();

    let runs: Collection<JobRunDoc> = get_cerebro_db_collection(&data, AdminCollection::Jobs);

    // Build a flexible filter: match by jid == {id} OR _id == {uuid}
    let filter = if let Ok(uuid) = Uuid::parse_str(&id) {
        doc! { "$or": [ { "id": uuid }, { "jid": &id } ] }
    } else {
        doc! { "jid": &id }
    };

    let doc_opt = match runs.find_one(filter).await {
        Ok(d) => d,
        Err(e) => {
            return HttpResponse::InternalServerError()
                .json(serde_json::json!({ "error": e.to_string() }))
        }
    };

    let Some(doc) = doc_opt else {
        return HttpResponse::NotFound().json(serde_json::json!({
            "completed": false,
            "status": "not_found"
        }));
    };

    // Normalize status and return a stable boolean + result/error
    let status = doc.status.clone();
    match status.as_str() {
        "succeeded" => HttpResponse::Ok().json(JobCompletionResponse {
            completed: true,
            result: doc.result,
            error: None,
            status: Some(status),
        }),
        "failed" => HttpResponse::Ok().json(JobCompletionResponse {
            completed: true,
            result: None,
            error: doc.error.or(Some("failed".into())),
            status: Some(status),
        }),
        "queued" | "running" | _ => HttpResponse::Ok().json(JobCompletionResponse {
            completed: false,
            result: None,
            error: None,
            status: Some(status),
        }),
    }
}

#[post("/jobs/schedule")]
pub async fn schedule_job_handler(
    data: web::Data<AppState>,
    body: web::Json<ScheduleJobRequest>,
    _: JwtAdminMiddleware,
) -> HttpResponse {
    let b = body.into_inner();

    let jobs_collection: Collection<ScheduleJob> =
        get_cerebro_db_collection(&data, AdminCollection::Scheduler);

    // Build strongly typed ScheduleJob
    let schedule = ScheduleJob::new(
        b.kind,
        b.args,
        b.queue,
        b.run_at,
        b.interval_seconds,
        b.retry,
        b.reserve_for_seconds,
        b.enabled,
    );
    let id = schedule.id.clone();

    match jobs_collection.insert_one(&schedule).await {
        Ok(_) => HttpResponse::Created().json(CreateScheduleResponse::success(&id)),
        Err(e) => HttpResponse::InternalServerError()
            .json(CreateScheduleResponse::server_error(e.to_string())),
    }
}


#[get("/jobs/schedule/status")]
pub async fn job_schedule_status_handler(
    data: web::Data<AppState>,
    _: JwtAdminMiddleware,
) -> HttpResponse {
    let jobs: Collection<ScheduleJob> = get_cerebro_db_collection(&data, AdminCollection::Scheduler);

    let now = chrono::Utc::now();

    // Next 10 jobs (enabled + upcoming)
    let mut next_cur = match jobs
        .find(doc! { "enabled": true, "run_at": { "$gte": now } })
        .sort(doc! { "run_at": 1i32 })
        .limit(10)
        .await
    {
        Ok(c) => c,
        Err(e) => return HttpResponse::InternalServerError()
            .json(JobsStatusResponse::server_error(e.to_string())),
    };

    let mut next = Vec::new();
    while let Some(job) = match next_cur.try_next().await {
        Ok(v) => v,
        Err(e) => return HttpResponse::InternalServerError()
            .json(JobsStatusResponse::server_error(e.to_string())),
    } {
        next.push(ScheduleJobSummary::from(&job));
    }

    // Previous 10 jobs (recently run)
    let mut prev_cur = match jobs
        .find(doc! { "last_run_at": { "$ne": Bson::Null } })
        .sort(doc! { "last_run_at": -1i32 })
        .limit(10)
        .await
    {
        Ok(c) => c,
        Err(e) => return HttpResponse::InternalServerError()
            .json(JobsStatusResponse::server_error(e.to_string())),
    };

    let mut previous = Vec::new();
    while let Some(job) = match prev_cur.try_next().await {
        Ok(v) => v,
        Err(e) => return HttpResponse::InternalServerError()
            .json(JobsStatusResponse::server_error(e.to_string())),
    } {
        previous.push(ScheduleJobSummary::from(&job));
    }

    if next.is_empty() && previous.is_empty() {
        return HttpResponse::Ok().json(JobsStatusResponse::not_found());
    }

    HttpResponse::Ok().json(JobsStatusResponse::success(JobsStatusData { next, previous }))
}


// Handler configuration
pub fn jobs_config(cfg: &mut web::ServiceConfig) {
    cfg.service(enqueue_now_handler)
        .service(schedule_job_handler)
        .service(job_completion_handler)
        .service(job_schedule_status_handler);
}