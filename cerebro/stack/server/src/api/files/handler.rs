use actix_web::{delete, get, patch, post, web, HttpResponse};
use cerebro_model::api::{
    files::{
        audit::{
            verify_chain, AuditAction, AuditActor, AuditError, AuditEvent, AUDIT_GENESIS_HASH,
        },
        response::{
            AuditTrailResponse, DeleteFileResponse, DeleteFilesResponse, ExpireResponse,
            GetFileResponse, ListFilesResponse, PurgeResponse, RegisterFileResponse,
            ReportOutResponse, RestoreTransitionResponse, UpdateLifecycleResponse,
            UpdateTagsResponse,
        },
        retention::{
            quarantine_grace_days_from_env, restore_available_days_from_env,
            warm_available_from_env, RestoreState, RetentionPolicy,
        },
        schema::{
            ExpireSchema, FileRelocateSchema, PurgeSchema, ReportOutSchema,
            RestoreTransitionSchema, UpdateFileLifecycleSchema, UpdateFileTagsSchema,
        },
    },
    teams::model::{Team, TeamAdminCollection},
};
use chrono::{DateTime, Duration, Utc};
use futures::TryStreamExt;
use mongodb::{
    bson::{doc, from_document, to_bson, Bson},
    error::{ErrorKind, WriteFailure},
    options::{FindOneOptions, IndexOptions},
    Collection, IndexModel,
};
use serde::Deserialize;
use std::collections::HashSet;
use std::sync::{Mutex, OnceLock};

use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::schema::RegisterFileSchema;
use cerebro_model::api::files::telemetry::{TelemetryEvent, TelemetryOp, TelemetryOutcome};

use crate::api::auth::jwt::{self, TeamAccessQuery};
use crate::api::server::AppState;

use crate::api::files::mongo::get_latest_files_paginated_pipeline;
use crate::api::utils::get_teams_db_collection;

#[post("/files/register")]
async fn register_file(
    data: web::Data<AppState>,
    schema: web::Json<RegisterFileSchema>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    match files_collection.find_one(doc! { "id": &schema.id }).await {
        Ok(None) => {}
        Ok(Some(_)) => {
            return HttpResponse::Conflict().json(RegisterFileResponse::conflict(&schema.id))
        }
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(RegisterFileResponse::server_error(err.to_string()))
        }
    }

    let result = files_collection
        .insert_one(SeaweedFile::from_schema(&schema))
        .await;

    match result {
        Ok(_) => {
            if let Err(err) = append_audit_event(
                &data,
                audit_team,
                AuditAction::Upload,
                Some(schema.id.clone()),
                schema.run_id.clone(),
                schema.sample_id.clone(),
                audit_actor,
                format!(
                    "uploaded '{}' ({} bytes, tier {:?}, retention {:?})",
                    schema.name, schema.size, schema.tier, schema.retention
                ),
            )
            .await
            {
                return HttpResponse::InternalServerError().json(
                    RegisterFileResponse::server_error(format!("audit append failed: {}", err)),
                );
            }
            data.metrics
                .record(&TelemetryEvent::success(TelemetryOp::Upload));
            HttpResponse::Ok().json(RegisterFileResponse::success(&schema.id))
        }
        Err(err) => HttpResponse::InternalServerError()
            .json(RegisterFileResponse::server_error(err.to_string())),
    }
}

#[derive(Deserialize)]
struct FileListQuery {
    // Paginated files
    page: u32,
    // Limit files returned - can be zero for all files
    limit: u32,
    // Optional run identifier
    run_id: Option<String>,
    // Optional biological sample identifier (linkage)
    sample_id: Option<String>,
    // Optional watcher identifier
    watcher_id: Option<String>,
}

#[get("/files")]
async fn list_files(
    data: web::Data<AppState>,
    query: web::Query<FileListQuery>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    let pipeline = get_latest_files_paginated_pipeline(
        query.run_id.clone(),
        query.sample_id.clone(),
        query.watcher_id.clone(),
        query.page as i64,
        query.limit as i64,
    );

    match files_collection.aggregate(pipeline).await {
        Ok(cursor) => {
            let files: Vec<SeaweedFile> = cursor
                .try_collect::<Vec<_>>()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .filter_map(|doc| from_document(doc).ok())
                .collect();

            match files.is_empty() {
                false => HttpResponse::Ok().json(ListFilesResponse::success(files)),
                true => HttpResponse::NotFound().json(ListFilesResponse::not_found()),
            }
        }
        Err(err) => HttpResponse::InternalServerError()
            .json(ListFilesResponse::server_error(err.to_string())),
    }
}

#[patch("/files/tags")]
async fn update_tags(
    data: web::Data<AppState>,
    schema: web::Json<UpdateFileTagsSchema>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    match files_collection
        .update_many(
            doc! { "id": { "$in": &schema.ids } },
            doc! { "$set": { "tags": to_bson(&schema.tags).unwrap() } },
        )
        .await
    {
        Ok(update_result) => {
            if update_result.matched_count > 0 {
                HttpResponse::Ok().json(UpdateTagsResponse::success())
            } else {
                HttpResponse::NotFound().json(UpdateTagsResponse::not_found())
            }
        }
        Err(err) => HttpResponse::InternalServerError()
            .json(UpdateTagsResponse::server_error(err.to_string())),
    }
}

#[derive(Deserialize)]
struct FilesDeleteQuery {
    // Optional run identifier
    run_id: Option<String>,
    // Optional sample identifier
    sample_id: Option<String>,
    // Delete all option
    all: Option<bool>,
}

#[delete("/files")]
async fn delete_files(
    data: web::Data<AppState>,
    query: web::Query<FilesDeleteQuery>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    let mut delete_query = doc! {}; // Default to deleting all if `all` is true

    if let Some(true) = query.all {
        // If `all` is true, delete all documents
        delete_query = doc! {};
    } else {
        if query.run_id.is_none() && query.sample_id.is_none() {
            return HttpResponse::BadRequest().json(DeleteFilesResponse::invalid_query());
        }

        // Build query based on provided parameters
        if let Some(run_id) = &query.run_id {
            delete_query.insert("run_id", run_id);
        }

        if let Some(sample_id) = &query.sample_id {
            delete_query.insert("sample_id", sample_id);
        }
    }

    // Protect held / still-retained files from bulk deletion: only delete
    // files that are not under legal hold and whose retention has lapsed (or was
    // never set). Protected files are left in place rather than failing the call.
    let now = Utc::now();
    delete_query.insert("legal_hold", doc! { "$ne": true });
    delete_query.insert(
        "$or",
        vec![
            Bson::Document(doc! { "retain_until": Bson::Null }),
            Bson::Document(doc! { "retain_until": { "$lte": to_bson(&now).unwrap() } }),
        ],
    );

    // Retrieve IDs of files to delete
    match files_collection.find(delete_query.clone()).await {
        Ok(cursor) => {
            let deleted_ids: Vec<String> = cursor
                .try_collect()
                .await
                .unwrap_or_else(|_| vec![])
                .into_iter()
                .map(|f| f.fid)
                .collect();

            // Proceed with deletion
            match files_collection.delete_many(delete_query).await {
                Ok(_) => {
                    let count = deleted_ids.len();
                    if let Err(err) = append_audit_event(
                        &data,
                        audit_team,
                        AuditAction::Delete,
                        None,
                        query.run_id.clone(),
                        query.sample_id.clone(),
                        audit_actor,
                        format!(
                            "bulk delete (all={:?}, run_id={:?}, sample_id={:?}, {} file(s))",
                            query.all, query.run_id, query.sample_id, count
                        ),
                    )
                    .await
                    {
                        return HttpResponse::InternalServerError().json(
                            DeleteFilesResponse::server_error(format!(
                                "audit append failed: {}",
                                err
                            )),
                        );
                    }
                    HttpResponse::Ok().json(DeleteFilesResponse::success(deleted_ids))
                }
                Err(err) => HttpResponse::InternalServerError()
                    .json(DeleteFilesResponse::server_error(err.to_string())),
            }
        }
        Err(err) => HttpResponse::InternalServerError()
            .json(DeleteFilesResponse::server_error(err.to_string())),
    }
}

#[delete("/files/{id}")]
async fn delete_file(
    data: web::Data<AppState>,
    id: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    let id = id.into_inner();

    // Enforce legal hold / active retention before deleting.
    match files_collection.find_one(doc! { "id": &id }).await {
        Ok(Some(file)) => {
            if let Some(reason) = protection_reason(&file, Utc::now()) {
                data.metrics
                    .record(&TelemetryEvent::rejected(TelemetryOp::Delete));
                return HttpResponse::Conflict().json(DeleteFileResponse::protected(&reason));
            }
        }
        Ok(None) => return HttpResponse::NotFound().json(DeleteFileResponse::not_found()),
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(DeleteFileResponse::server_error(err.to_string()))
        }
    }

    match files_collection
        .find_one_and_delete(doc! { "id":  &id})
        .await
    {
        Ok(deleted) => match deleted {
            Some(file) => {
                if let Err(err) = append_audit_event(
                    &data,
                    audit_team,
                    AuditAction::Delete,
                    Some(file.id.clone()),
                    file.run_id.clone(),
                    file.sample_id.clone(),
                    audit_actor,
                    format!("deleted '{}' (fid {})", file.name, file.fid),
                )
                .await
                {
                    return HttpResponse::InternalServerError().json(
                        DeleteFileResponse::server_error(format!("audit append failed: {}", err)),
                    );
                }
                data.metrics
                    .record(&TelemetryEvent::success(TelemetryOp::Delete));
                HttpResponse::Ok().json(DeleteFileResponse::success(file))
            }
            None => HttpResponse::NotFound().json(DeleteFileResponse::not_found()),
        },
        Err(err) => HttpResponse::InternalServerError()
            .json(DeleteFileResponse::server_error(err.to_string())),
    }
}

// Handler configuration
/// Persist a partial lifecycle update to a single file.
///
/// Backs the report-out (FS-7/FS-9), restore (FS-4), and integrity (FS-6) flows:
/// each supplies only the fields it computed. Empty updates are rejected; only
/// the supplied fields are `$set`, leaving the rest (and concurrent tag updates)
/// untouched.
#[patch("/files/{id}/lifecycle")]
async fn update_file_lifecycle(
    data: web::Data<AppState>,
    id: web::Path<String>,
    schema: web::Json<UpdateFileLifecycleSchema>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let id = id.into_inner();

    if schema.is_empty() {
        return HttpResponse::BadRequest().json(UpdateLifecycleResponse::invalid_query());
    }

    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    // Build a `$set` document from only the supplied fields.
    let mut set = doc! {};
    if let Some(tier) = &schema.tier {
        set.insert("tier", to_bson(tier).unwrap());
        // Committing a tier stamps the move time and clears any in-flight claim.
        set.insert("tier_moved_at", to_bson(&Utc::now()).unwrap());
        set.insert("pending_tier", Bson::Null);
        set.insert("pending_since", Bson::Null);
    }
    if schema.clear_retain_until {
        set.insert("retain_until", Bson::Null);
    } else if let Some(retain_until) = &schema.retain_until {
        set.insert("retain_until", to_bson(retain_until).unwrap());
    }
    if let Some(reported_at) = &schema.reported_at {
        set.insert("reported_at", to_bson(reported_at).unwrap());
    }
    if let Some(legal_hold) = schema.legal_hold {
        set.insert("legal_hold", legal_hold);
    }
    if let Some(archived) = schema.archived {
        set.insert("archived", archived);
    }
    if let Some(restore_state) = &schema.restore_state {
        set.insert("restore_state", to_bson(restore_state).unwrap());
    }
    if let Some(replicas) = schema.replicas {
        set.insert("replicas", replicas as i64);
    }
    // Pending-tier claim (only when not committing a tier, which auto-clears it).
    if schema.tier.is_none() {
        if schema.clear_pending_tier {
            set.insert("pending_tier", Bson::Null);
        } else if let Some(pending_tier) = &schema.pending_tier {
            set.insert("pending_tier", to_bson(pending_tier).unwrap());
        }
        // Claim timestamp: stamped with the claim, cleared on abort.
        if schema.clear_pending_since {
            set.insert("pending_since", Bson::Null);
        } else if let Some(pending_since) = &schema.pending_since {
            set.insert("pending_since", to_bson(pending_since).unwrap());
        }
    }
    // Verification timestamp: stamped by a successful verify.
    if let Some(verified_at) = &schema.verified_at {
        set.insert("verified_at", to_bson(verified_at).unwrap());
    }

    // Classify the change for the audit trail and summarise the supplied fields.
    let action = if schema.legal_hold.is_some() {
        AuditAction::LegalHoldChange
    } else if schema.restore_state.is_some() {
        AuditAction::Restore
    } else if schema.archived.is_some() || schema.tier.is_some() {
        AuditAction::TierMove
    } else if schema.reported_at.is_some() {
        AuditAction::ReportOut
    } else {
        AuditAction::TierMove
    };
    let mut changes: Vec<String> = Vec::new();
    if let Some(tier) = &schema.tier {
        changes.push(format!("tier={:?}", tier));
    }
    if schema.clear_retain_until {
        changes.push("retain_until=cleared".into());
    } else if let Some(ru) = &schema.retain_until {
        changes.push(format!("retain_until={}", ru.to_rfc3339()));
    }
    if let Some(ra) = &schema.reported_at {
        changes.push(format!("reported_at={}", ra.to_rfc3339()));
    }
    if let Some(lh) = schema.legal_hold {
        changes.push(format!("legal_hold={}", lh));
    }
    if let Some(ar) = schema.archived {
        changes.push(format!("archived={}", ar));
    }
    if let Some(rs) = &schema.restore_state {
        changes.push(format!("restore_state={:?}", rs));
    }
    if let Some(rp) = schema.replicas {
        changes.push(format!("replicas={}", rp));
    }
    if schema.clear_pending_tier {
        changes.push("pending_tier=cleared".into());
    } else if let Some(pt) = &schema.pending_tier {
        changes.push(format!("pending_tier={:?}", pt));
    }
    if schema.clear_pending_since {
        changes.push("pending_since=cleared".into());
    } else if let Some(ps) = &schema.pending_since {
        changes.push(format!("pending_since={}", ps.to_rfc3339()));
    }
    if let Some(va) = &schema.verified_at {
        changes.push(format!("verified_at={}", va.to_rfc3339()));
    }
    if let Some(et) = &schema.expected_tier {
        changes.push(format!("expected_tier={:?}", et));
    }
    let detail = changes.join(", ");

    // Telemetry: capture op + bounded tier detail before `action`/`detail`
    // are moved into the audit append below.
    let tel_op = telemetry_op_for(&action);
    let tel_detail = schema
        .tier
        .as_ref()
        .map(|t| format!("{t:?}").to_lowercase());

    // Compare-and-set: when expected_tier is supplied, the update applies only if
    // the file's current tier matches — making mover claims/commits idempotent.
    let mut filter = doc! { "id": &id };
    if let Some(expected) = &schema.expected_tier {
        filter.insert("tier", to_bson(expected).unwrap());
    }

    match files_collection
        .update_one(filter, doc! { "$set": set })
        .await
    {
        Ok(update_result) => {
            if update_result.matched_count > 0 {
                if let Err(err) = append_audit_event(
                    &data,
                    audit_team,
                    action,
                    Some(id.clone()),
                    None,
                    None,
                    audit_actor,
                    detail,
                )
                .await
                {
                    return HttpResponse::InternalServerError().json(
                        UpdateLifecycleResponse::server_error(format!(
                            "audit append failed: {}",
                            err
                        )),
                    );
                }
                if let Some(op) = tel_op {
                    let ev = match tel_detail {
                        Some(d) => TelemetryEvent::with_detail(op, TelemetryOutcome::Success, d),
                        None => TelemetryEvent::success(op),
                    };
                    data.metrics.record(&ev);
                }
                HttpResponse::Ok().json(UpdateLifecycleResponse::success(&id))
            } else if schema.expected_tier.is_some() {
                // No match with a precondition set: distinguish a missing file from
                // an unmet tier precondition (an idempotent no-op for a mover).
                match files_collection.find_one(doc! { "id": &id }).await {
                    Ok(Some(_)) => {
                        // Compare-and-set rejected: a benign, expected mover no-op.
                        data.metrics.record(&TelemetryEvent::rejected(
                            tel_op.unwrap_or(TelemetryOp::TierMove),
                        ));
                        HttpResponse::Conflict()
                            .json(UpdateLifecycleResponse::precondition_failed(&id))
                    }
                    Ok(None) => {
                        HttpResponse::NotFound().json(UpdateLifecycleResponse::not_found(&id))
                    }
                    Err(err) => HttpResponse::InternalServerError()
                        .json(UpdateLifecycleResponse::server_error(err.to_string())),
                }
            } else {
                HttpResponse::NotFound().json(UpdateLifecycleResponse::not_found(&id))
            }
        }
        Err(err) => HttpResponse::InternalServerError()
            .json(UpdateLifecycleResponse::server_error(err.to_string())),
    }
}

/// Fetch a single file record by id.
///
/// The frontend uses this to resolve an artefact's metadata and lifecycle state
/// (tier, retain_until, reported_at, legal_hold, restore_state, replicas) before
/// initiating a download.
#[get("/files/{id}")]
async fn get_file(
    data: web::Data<AppState>,
    id: web::Path<String>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let id = id.into_inner();
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    match files_collection.find_one(doc! { "id": &id }).await {
        Ok(Some(file)) => HttpResponse::Ok().json(GetFileResponse::success(file)),
        Ok(None) => HttpResponse::NotFound().json(GetFileResponse::not_found(&id)),
        Err(err) => {
            HttpResponse::InternalServerError().json(GetFileResponse::server_error(err.to_string()))
        }
    }
}

/// Report-out trigger: anchor retention and move a run's files off the hot tier.
///
/// The authoritative entry point fired when a result is reported out (by report
/// generation or an operator). The server applies its **configured** retention
/// policy ([`RetentionPolicy::from_env`]) and warm-tier availability
/// (`CEREBRO_FS_WARM_AVAILABLE`) to every matching file, persisting the
/// re-anchored `retain_until`, the `reported_at` anchor, and the new tier
/// (warm when available, else cold). The physical hot→warm→cold(S3) move is the
/// worker's job; this records the decision.
#[post("/files/report-out")]
async fn report_out_files(
    data: web::Data<AppState>,
    schema: web::Json<ReportOutSchema>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    let reported_at = schema.reported_at.unwrap_or_else(Utc::now);
    let policy = RetentionPolicy::from_env();
    let warm_available = warm_available_from_env();

    let mut filter = doc! { "run_id": &schema.run_id };
    if let Some(sample_id) = &schema.sample_id {
        filter.insert("sample_id", sample_id);
    }

    let mut cursor = match files_collection.find(filter).await {
        Ok(cursor) => cursor,
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(ReportOutResponse::server_error(err.to_string()))
        }
    };

    let mut updated: u64 = 0;
    loop {
        match cursor.try_next().await {
            Ok(Some(file)) => {
                let transition =
                    policy.report_out_transition(file.retention, reported_at, warm_available);
                let audit_detail = format!(
                    "reported out -> tier {:?}, retain_until {:?}",
                    transition.target_tier, transition.retain_until
                );

                let mut set = doc! {
                    "tier": to_bson(&transition.target_tier).unwrap(),
                    "reported_at": to_bson(&transition.reported_at).unwrap(),
                };
                match transition.retain_until {
                    Some(retain_until) => {
                        set.insert("retain_until", to_bson(&retain_until).unwrap());
                    }
                    None => {
                        set.insert("retain_until", Bson::Null);
                    }
                }

                if let Err(err) = files_collection
                    .update_one(doc! { "id": &file.id }, doc! { "$set": set })
                    .await
                {
                    return HttpResponse::InternalServerError()
                        .json(ReportOutResponse::server_error(err.to_string()));
                }
                if let Err(err) = append_audit_event(
                    &data,
                    audit_team.clone(),
                    AuditAction::ReportOut,
                    Some(file.id.clone()),
                    file.run_id.clone(),
                    file.sample_id.clone(),
                    audit_actor.clone(),
                    audit_detail,
                )
                .await
                {
                    return HttpResponse::InternalServerError().json(
                        ReportOutResponse::server_error(format!("audit append failed: {}", err)),
                    );
                }
                data.metrics.record(&TelemetryEvent::with_detail(
                    TelemetryOp::ReportOut,
                    TelemetryOutcome::Success,
                    format!("{:?}", transition.target_tier).to_lowercase(),
                ));
                updated += 1;
            }
            Ok(None) => break,
            Err(err) => {
                return HttpResponse::InternalServerError()
                    .json(ReportOutResponse::server_error(err.to_string()))
            }
        }
    }

    if updated == 0 {
        return HttpResponse::NotFound().json(ReportOutResponse::not_found(&schema.run_id));
    }
    HttpResponse::Ok().json(ReportOutResponse::success(
        updated,
        reported_at.to_rfc3339(),
    ))
}

/// Build an [`AuditActor`] from the authenticated principal.
fn actor_of(auth_guard: &jwt::JwtDataMiddleware) -> AuditActor {
    AuditActor {
        id: auth_guard.user.id.clone(),
        email: auth_guard.user.email.clone(),
    }
}

/// Map an [`AuditAction`] to a telemetry op. `None` for actions we don't
/// surface as lifecycle-movement metrics (legal-hold / tag edits).
fn telemetry_op_for(action: &AuditAction) -> Option<TelemetryOp> {
    match action {
        AuditAction::Upload => Some(TelemetryOp::Upload),
        AuditAction::TierMove => Some(TelemetryOp::TierMove),
        AuditAction::ReportOut => Some(TelemetryOp::ReportOut),
        AuditAction::Restore => Some(TelemetryOp::Restore),
        AuditAction::Expiry => Some(TelemetryOp::Expire),
        AuditAction::Delete => Some(TelemetryOp::Delete),
        AuditAction::LegalHoldChange | AuditAction::TagChange => None,
    }
}

/// Append a sealed event to the team's audit chain.
///
/// Reads the current chain tip (highest sequence), links the new event to it,
/// seals it (BLAKE3 over `prev_hash` + body) and inserts it. Best-effort: a
/// failure is logged but does not fail the triggering operation (see the
/// notes on fail-open vs fail-closed auditing).
#[allow(clippy::too_many_arguments)]
/// Reason a file must not be deleted/purged, or `None` when deletion is allowed.
///
/// Enforces the two compliance invariants server-side: a file under legal
/// hold, or still within its retention period, is protected. Clearing the hold or
/// retention via the lifecycle endpoint (an audited admin action) is the
/// sanctioned way to release it.
fn protection_reason(file: &SeaweedFile, now: DateTime<Utc>) -> Option<String> {
    if file.legal_hold {
        return Some(format!("file {} is under legal hold", file.id));
    }
    if let Some(until) = file.retain_until {
        if until > now {
            return Some(format!(
                "file {} is retained until {}",
                file.id,
                until.to_rfc3339()
            ));
        }
    }
    None
}

/// Maximum re-reads of the chain tip when concurrent appends contend for the
/// same sequence number.
const AUDIT_MAX_RETRIES: u32 = 8;

/// Whether audit appends are best-effort (`CEREBRO_AUDIT_FAIL_OPEN`). Default is
/// fail-closed: an unauditable lifecycle action surfaces as an operation failure.
fn audit_fail_open() -> bool {
    matches!(
        std::env::var("CEREBRO_AUDIT_FAIL_OPEN").ok().as_deref(),
        Some("true") | Some("1")
    )
}

/// Per-process set of team chains whose unique `sequence` index has been ensured.
fn audit_index_ensured() -> &'static Mutex<HashSet<String>> {
    static ENSURED: OnceLock<Mutex<HashSet<String>>> = OnceLock::new();
    ENSURED.get_or_init(|| Mutex::new(HashSet::new()))
}

/// True for a MongoDB duplicate-key error (E11000) — a lost sequence race.
fn is_duplicate_key(err: &mongodb::error::Error) -> bool {
    matches!(&*err.kind, ErrorKind::Write(WriteFailure::WriteError(we)) if we.code == 11000)
}

/// Ensure the unique index on `sequence` for a team's audit chain (idempotent;
/// ensured once per team per process). The index is what serialises concurrent
/// appends: a lost race fails with a duplicate key and is retried.
async fn ensure_audit_index(collection: &Collection<AuditEvent>, team_key: &str) {
    if audit_index_ensured().lock().unwrap().contains(team_key) {
        return;
    }
    let model = IndexModel::builder()
        .keys(doc! { "sequence": 1 })
        .options(IndexOptions::builder().unique(true).build())
        .build();
    match collection.create_index(model).await {
        Ok(_) => {
            audit_index_ensured()
                .lock()
                .unwrap()
                .insert(team_key.to_string());
        }
        Err(err) => log::error!("Failed to ensure audit sequence index: {}", err),
    }
}

/// Apply the fail-open/closed policy to a terminal append failure.
fn audit_failure(err: AuditError) -> Result<(), AuditError> {
    if audit_fail_open() {
        log::error!("Audit append failed (fail-open): {}", err);
        Ok(())
    } else {
        Err(err)
    }
}

/// Append a sealed event to the team's audit chain.
///
/// Ensures the unique `sequence` index, then reads the chain tip, links + seals
/// the new event and inserts it. On a duplicate-key (lost sequence race) the tip
/// is re-read and the append retried up to [`AUDIT_MAX_RETRIES`], so concurrent
/// appends serialise rather than fork. Returns `Err` on terminal failure unless
/// `CEREBRO_AUDIT_FAIL_OPEN` relaxes to best-effort.
#[allow(clippy::too_many_arguments)]
async fn append_audit_event(
    data: &web::Data<AppState>,
    team: Team,
    action: AuditAction,
    file_id: Option<String>,
    run_id: Option<String>,
    sample_id: Option<String>,
    actor: AuditActor,
    detail: String,
) -> Result<(), AuditError> {
    let team_key = team.id.clone();
    let audit_collection: Collection<AuditEvent> =
        get_teams_db_collection(data, team, TeamAdminCollection::AuditLogs);

    ensure_audit_index(&audit_collection, &team_key).await;

    let mut attempt: u32 = 0;
    loop {
        attempt += 1;

        let (sequence, prev_hash) = match audit_collection
            .find_one(doc! {})
            .with_options(
                FindOneOptions::builder()
                    .sort(doc! { "sequence": -1 })
                    .build(),
            )
            .await
        {
            Ok(Some(tip)) => (tip.sequence + 1, tip.hash),
            Ok(None) => (0, AUDIT_GENESIS_HASH.to_string()),
            Err(err) => {
                data.metrics
                    .record(&TelemetryEvent::failure(TelemetryOp::Audit));
                return audit_failure(AuditError::TipLookup(err.to_string()));
            }
        };

        let mut event = AuditEvent::new(
            sequence,
            Utc::now(),
            action.clone(),
            file_id.clone(),
            run_id.clone(),
            sample_id.clone(),
            actor.clone(),
            detail.clone(),
            prev_hash,
        );
        event.seal();

        match audit_collection.insert_one(event).await {
            Ok(_) => {
                data.metrics
                    .record(&TelemetryEvent::success(TelemetryOp::Audit));
                return Ok(());
            }
            Err(err) if is_duplicate_key(&err) => {
                if attempt >= AUDIT_MAX_RETRIES {
                    data.metrics
                        .record(&TelemetryEvent::failure(TelemetryOp::Audit));
                    return audit_failure(AuditError::RetriesExhausted);
                }
                // Lost the sequence race; re-read the tip and try again.
                continue;
            }
            Err(err) => {
                data.metrics
                    .record(&TelemetryEvent::failure(TelemetryOp::Audit));
                return audit_failure(AuditError::Insert(err.to_string()));
            }
        }
    }
}

#[derive(Deserialize)]
struct AuditQuery {
    // Optional file identifier to filter the returned events
    file_id: Option<String>,
    // Optional run identifier to filter the returned events
    run_id: Option<String>,
}

/// Return the audit trail (optionally filtered) plus whole-chain verification.
///
/// The full team chain is fetched in sequence order and verified for integrity
/// (`verified`); the returned `events` are then narrowed by any `file_id`/`run_id`
/// filter. Verification is always over the entire chain, so a filtered view still
/// reports whether the underlying trail has been tampered with.
#[get("/audit")]
async fn get_audit(
    data: web::Data<AppState>,
    query: web::Query<AuditQuery>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let audit_collection: Collection<AuditEvent> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::AuditLogs);

    let cursor = match audit_collection
        .find(doc! {})
        .sort(doc! { "sequence": 1 })
        .await
    {
        Ok(cursor) => cursor,
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(AuditTrailResponse::server_error(err.to_string()))
        }
    };
    let all: Vec<AuditEvent> = cursor.try_collect().await.unwrap_or_else(|_| vec![]);

    let verified = verify_chain(&all);

    let events: Vec<AuditEvent> = all
        .into_iter()
        .filter(|e| match &query.file_id {
            Some(fid) => e.file_id.as_deref() == Some(fid.as_str()),
            None => true,
        })
        .filter(|e| match &query.run_id {
            Some(rid) => e.run_id.as_deref() == Some(rid.as_str()),
            None => true,
        })
        .collect();

    if events.is_empty() {
        return HttpResponse::NotFound().json(AuditTrailResponse::not_found());
    }
    HttpResponse::Ok().json(AuditTrailResponse::success(events, verified))
}

/// Non-destructive expiry sweep: quarantine files past retention.
///
/// Selects files whose `retain_until` has lapsed, that are **not** under legal
/// hold, and that are still active, optionally scoped to a run/sample. Each is
/// moved to `ExpiryState::Quarantined` (with `quarantined_at`) and audited — never
/// deleted. `dry_run` reports the eligible count without changing state. The
/// hard removal of quarantined data is the separate, gated `purge`.
#[post("/files/expire")]
async fn expire_files(
    data: web::Data<AppState>,
    schema: web::Json<ExpireSchema>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    let now = Utc::now();
    let mut filter = doc! {
        "legal_hold": { "$ne": true },
        "retain_until": { "$lte": to_bson(&now).unwrap(), "$ne": Bson::Null },
        "expiry_state": { "$ne": "quarantined" },
    };
    if let Some(run_id) = &schema.run_id {
        filter.insert("run_id", run_id);
    }
    if let Some(sample_id) = &schema.sample_id {
        filter.insert("sample_id", sample_id);
    }

    let mut cursor = match files_collection.find(filter).await {
        Ok(cursor) => cursor,
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(ExpireResponse::server_error(err.to_string()))
        }
    };

    let mut quarantined: u64 = 0;
    loop {
        match cursor.try_next().await {
            Ok(Some(file)) => {
                if !schema.dry_run {
                    let set = doc! { "expiry_state": "quarantined", "quarantined_at": to_bson(&now).unwrap() };
                    if let Err(err) = files_collection
                        .update_one(doc! { "id": &file.id }, doc! { "$set": set })
                        .await
                    {
                        return HttpResponse::InternalServerError()
                            .json(ExpireResponse::server_error(err.to_string()));
                    }
                    if let Err(err) = append_audit_event(
                        &data,
                        audit_team.clone(),
                        AuditAction::Expiry,
                        Some(file.id.clone()),
                        file.run_id.clone(),
                        file.sample_id.clone(),
                        audit_actor.clone(),
                        format!("quarantined (retain_until {:?} lapsed)", file.retain_until),
                    )
                    .await
                    {
                        return HttpResponse::InternalServerError().json(
                            ExpireResponse::server_error(format!("audit append failed: {}", err)),
                        );
                    }
                    data.metrics
                        .record(&TelemetryEvent::success(TelemetryOp::Expire));
                }
                quarantined += 1;
            }
            Ok(None) => break,
            Err(err) => {
                return HttpResponse::InternalServerError()
                    .json(ExpireResponse::server_error(err.to_string()))
            }
        }
    }

    HttpResponse::Ok().json(ExpireResponse::success(quarantined, schema.dry_run))
}

/// Gated hard purge: permanently remove quarantined files past the grace window
///.
///
/// Selects files in `ExpiryState::Quarantined`, **not** under legal hold, whose
/// `quarantined_at` is older than the configured grace window
/// (`CEREBRO_QUARANTINE_GRACE_DAYS`). Each record is deleted and audited; the
/// returned `fids` let a worker reclaim the bytes. `dry_run` reports the eligible
/// set without deleting anything.
#[post("/files/purge")]
async fn purge_files(
    data: web::Data<AppState>,
    schema: web::Json<PurgeSchema>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    let now = Utc::now();
    let grace = quarantine_grace_days_from_env();
    let cutoff = now - Duration::days(grace);

    let mut filter = doc! {
        "legal_hold": { "$ne": true },
        "expiry_state": "quarantined",
        "quarantined_at": { "$lte": to_bson(&cutoff).unwrap(), "$ne": Bson::Null },
    };
    if let Some(run_id) = &schema.run_id {
        filter.insert("run_id", run_id);
    }
    if let Some(sample_id) = &schema.sample_id {
        filter.insert("sample_id", sample_id);
    }

    let mut cursor = match files_collection.find(filter).await {
        Ok(cursor) => cursor,
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(PurgeResponse::server_error(err.to_string()))
        }
    };

    let mut fids: Vec<String> = Vec::new();
    loop {
        match cursor.try_next().await {
            Ok(Some(file)) => {
                if !schema.dry_run {
                    if let Err(err) = files_collection.delete_one(doc! { "id": &file.id }).await {
                        return HttpResponse::InternalServerError()
                            .json(PurgeResponse::server_error(err.to_string()));
                    }
                    if let Err(err) = append_audit_event(
                        &data,
                        audit_team.clone(),
                        AuditAction::Delete,
                        Some(file.id.clone()),
                        file.run_id.clone(),
                        file.sample_id.clone(),
                        audit_actor.clone(),
                        format!("purged after quarantine grace (fid {})", file.fid),
                    )
                    .await
                    {
                        return HttpResponse::InternalServerError().json(
                            PurgeResponse::server_error(format!("audit append failed: {}", err)),
                        );
                    }
                    data.metrics
                        .record(&TelemetryEvent::success(TelemetryOp::Purge));
                }
                fids.push(file.fid.clone());
            }
            Ok(None) => break,
            Err(err) => {
                return HttpResponse::InternalServerError()
                    .json(PurgeResponse::server_error(err.to_string()))
            }
        }
    }

    HttpResponse::Ok().json(PurgeResponse::success(fids, schema.dry_run))
}

/// Drive a restore state-machine transition.
///
/// Validates the requested `target` against the current `restore_state` via the
/// model's transition guard, applies it with a compare-and-set on the observed
/// state (atomic against concurrent workers), and stamps the relevant timestamps
/// (`restore_requested_at` on `Requested`; `restore_available_at` +
/// `restore_expires_at` = now + `CEREBRO_RESTORE_AVAILABLE_DAYS` on `Restored`).
/// A no-op (already in `target`) returns success without resetting timestamps; an
/// invalid transition returns `409`. The restore worker drives this on poll.
#[post("/files/{id}/restore")]
async fn restore_transition(
    data: web::Data<AppState>,
    id: web::Path<String>,
    schema: web::Json<RestoreTransitionSchema>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let id = id.into_inner();
    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    let target = schema.target;

    let file = match files_collection.find_one(doc! { "id": &id }).await {
        Ok(Some(file)) => file,
        Ok(None) => {
            return HttpResponse::NotFound().json(RestoreTransitionResponse::not_found(&id))
        }
        Err(err) => {
            return HttpResponse::InternalServerError()
                .json(RestoreTransitionResponse::server_error(err.to_string()))
        }
    };
    let current = file.restore_state;

    // Idempotent no-op: already in the target state (don't reset timestamps).
    if current == target {
        return HttpResponse::Ok().json(RestoreTransitionResponse::noop(current));
    }
    if !current.can_transition_to(&target) {
        data.metrics
            .record(&TelemetryEvent::rejected(TelemetryOp::Restore));
        return HttpResponse::Conflict().json(RestoreTransitionResponse::invalid_transition(
            &current.to_string(),
            &target.to_string(),
        ));
    }

    let now = Utc::now();
    let mut set = doc! { "restore_state": to_bson(&target).unwrap() };
    let mut expires_at = None;
    match &target {
        RestoreState::Requested => {
            set.insert("restore_requested_at", to_bson(&now).unwrap());
        }
        RestoreState::Restored => {
            let exp = now + Duration::days(restore_available_days_from_env());
            set.insert("restore_available_at", to_bson(&now).unwrap());
            set.insert("restore_expires_at", to_bson(&exp).unwrap());
            expires_at = Some(exp);
        }
        _ => {}
    }

    // Compare-and-set on the observed state: rejects a concurrent advance.
    let filter = doc! { "id": &id, "restore_state": to_bson(&current).unwrap() };
    match files_collection
        .update_one(filter, doc! { "$set": set })
        .await
    {
        Ok(result) if result.matched_count > 0 => {
            if let Err(err) = append_audit_event(
                &data,
                audit_team,
                AuditAction::Restore,
                Some(id.clone()),
                file.run_id.clone(),
                file.sample_id.clone(),
                audit_actor,
                format!("restore {} -> {}", current, target),
            )
            .await
            {
                return HttpResponse::InternalServerError().json(
                    RestoreTransitionResponse::server_error(format!(
                        "audit append failed: {}",
                        err
                    )),
                );
            }
            // Kick off the restore worker promptly when a restore is (re)requested
            //. Best-effort: the hourly `restore_scan` re-drives anything
            // that slips, so a transient Faktory hiccup here is not fatal.
            if target == RestoreState::Requested {
                match faktory::Client::connect().await {
                    Ok(mut fk) => {
                        let job = faktory::Job::new(
                            "restore_drive",
                            vec![serde_json::json!({ "file_id": id.clone() })],
                        )
                        .on_queue("lifecycle");
                        if let Err(e) = fk.enqueue(job).await {
                            log::warn!("restore requested for {id} but failed to enqueue restore_drive (restore_scan will retry): {e}");
                        }
                    }
                    Err(e) => log::warn!("restore requested for {id} but Faktory unavailable to enqueue restore_drive (restore_scan will retry): {e}"),
                }
            }
            data.metrics.record(&TelemetryEvent::with_detail(
                TelemetryOp::Restore,
                TelemetryOutcome::Success,
                target.to_string().to_lowercase(),
            ));
            HttpResponse::Ok().json(RestoreTransitionResponse::success(target, expires_at))
        }
        Ok(_) => {
            data.metrics
                .record(&TelemetryEvent::rejected(TelemetryOp::Restore));
            HttpResponse::Conflict().json(RestoreTransitionResponse::invalid_transition(
                &current.to_string(),
                &target.to_string(),
            ))
        }
        Err(err) => HttpResponse::InternalServerError()
            .json(RestoreTransitionResponse::server_error(err.to_string())),
    }
}

/// Dedicated archive/relocate endpoint.
///
/// Repoints a file between local and remote (archival) storage in one
/// compare-and-set: it sets `tier`, `archived`, optionally a new `fid`, and the
/// `archive_key`, applied only if the file's current `tier` equals
/// `expected_tier`. Kept separate from the lifecycle update so a fid repoint is
/// never entangled with routine edits. Records a `TierMove` audit event with an
/// `archive`/`restore` detail.
#[post("/files/{id}/relocate")]
async fn relocate_file(
    data: web::Data<AppState>,
    id: web::Path<String>,
    schema: web::Json<FileRelocateSchema>,
    _: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware,
) -> HttpResponse {
    let id = id.into_inner();
    let audit_team = auth_guard.team.clone();
    let audit_actor = actor_of(&auth_guard);
    let files_collection: Collection<SeaweedFile> =
        get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Files);

    // A relocation always commits a tier, so it stamps tier_moved_at and clears
    // any in-flight claim — same invariant the lifecycle tier-commit upholds.
    let mut set = doc! {
        "tier": to_bson(&schema.target_tier).unwrap(),
        "tier_moved_at": to_bson(&Utc::now()).unwrap(),
        "pending_tier": Bson::Null,
        "pending_since": Bson::Null,
        "archived": schema.archived,
    };
    if let Some(fid) = &schema.fid {
        set.insert("fid", fid.clone());
    }
    if schema.clear_archive_key {
        set.insert("archive_key", Bson::Null);
    } else if let Some(key) = &schema.archive_key {
        set.insert("archive_key", key.clone());
    }

    let detail = if schema.archived {
        "archive"
    } else {
        "restore"
    };

    // CAS on the current tier so concurrent movers can't double-apply.
    let filter = doc! { "id": &id, "tier": to_bson(&schema.expected_tier).unwrap() };

    match files_collection
        .update_one(filter, doc! { "$set": set })
        .await
    {
        Ok(result) => {
            if result.matched_count > 0 {
                if let Err(err) = append_audit_event(
                    &data,
                    audit_team,
                    AuditAction::TierMove,
                    Some(id.clone()),
                    None,
                    None,
                    audit_actor,
                    detail.to_string(),
                )
                .await
                {
                    return HttpResponse::InternalServerError().json(
                        UpdateLifecycleResponse::server_error(format!(
                            "audit append failed: {}",
                            err
                        )),
                    );
                }
                data.metrics.record(&TelemetryEvent::with_detail(
                    TelemetryOp::TierMove,
                    TelemetryOutcome::Success,
                    detail.to_string(),
                ));
                HttpResponse::Ok().json(UpdateLifecycleResponse::success(&id))
            } else {
                // Distinguish a missing file from an unmet CAS precondition.
                match files_collection.find_one(doc! { "id": &id }).await {
                    Ok(Some(_)) => {
                        data.metrics
                            .record(&TelemetryEvent::rejected(TelemetryOp::TierMove));
                        HttpResponse::Conflict()
                            .json(UpdateLifecycleResponse::precondition_failed(&id))
                    }
                    Ok(None) => {
                        HttpResponse::NotFound().json(UpdateLifecycleResponse::not_found(&id))
                    }
                    Err(err) => HttpResponse::InternalServerError()
                        .json(UpdateLifecycleResponse::server_error(err.to_string())),
                }
            }
        }
        Err(err) => HttpResponse::InternalServerError()
            .json(UpdateLifecycleResponse::server_error(err.to_string())),
    }
}

pub fn files_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_file)
        .service(delete_file)
        .service(delete_files)
        .service(update_tags)
        .service(update_file_lifecycle)
        .service(relocate_file)
        .service(report_out_files)
        .service(expire_files)
        .service(purge_files)
        .service(restore_transition)
        .service(get_audit)
        .service(get_file)
        .service(list_files);
}
