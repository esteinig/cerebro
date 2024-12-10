use actix_web::body::MessageBody;
use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;
use futures::stream::TryStreamExt;
use mongodb::{bson::doc, Collection};
use base64::engine::general_purpose;
use base64::Engine;
use actix_web::{get, post, web, HttpResponse, patch, delete, HttpRequest};

use crate::api::auth::jwt;
use cerebro_model::api::config::Config;
use cerebro_report::report::{TemplateConfig, AssayTemplate, BioinformaticsTemplate};
use cerebro_model::api::users::model::Role;
use crate::api::utils::{as_csv_string, get_teams_db_collection};
use crate::api::server::AppState;
use crate::api::cerebro::mongo::*;
use crate::api::logs::utils::log_database_change;
use cerebro_model::api::logs::model::{LogModule, Action, RequestLog, AccessDetails};
use cerebro_model::api::teams::model::{DatabaseId, ProjectCollection, ProjectId, Team, TeamAdminCollection, TeamDatabase};
use cerebro_model::api::cerebro::model::{
    Cerebro, 
    PriorityTaxon, 
    ReportEntry, 
    SampleComment, 
    WorkflowId, 
    CerebroId, 
    DecisionType, 
    PriorityTaxonDecision, 
    SampleId,
    SampleConfig,
    WorkflowConfig,
    RunConfig
};
use cerebro_model::api::cerebro::schema::{
    PriorityTaxonDecisionSchema,
    SampleCommentSchema,
    SampleDeleteSchema,
    SampleSummaryQcSchema,
    SampleDescriptionSchema,
    TaxaSummarySchema,
    PriorityTaxonSchema,
    ReportSchema
};
use cerebro_model::api::cerebro::response::{TaxonSummaryOverview, TaxaSummaryMongoPipeline};
use cerebro_report::report::ClinicalReport;

use cerebro_pipe::taxa::filters::*;
use cerebro_pipe::modules::quality::{QualityControl, ReadQualityControl, ModelConfig};
use cerebro_pipe::taxa::taxon::{Taxon, TaxonOverview, aggregate};

type CerebroIds = String;

fn qc_config_from_model(sample_config: Option<SampleConfig>, run_config: Option<RunConfig>, workflow_config: Option<WorkflowConfig>) -> ModelConfig {

    let (sample_type, sample_group, ercc_input_mass, library_tag) = match sample_config {
        Some(config) => (Some(config.sample_type), Some(config.sample_group), config.ercc_input_mass, Some(config.tags.join("-"))), 
        None => (None, None, None, None)
    };

    let (run_id, run_date) = match run_config {
        Some(config) => (Some(config.id), Some(config.date)),
        None => (None, None)
    };

    let (workflow_name, workflow_id, workflow_date) = match workflow_config {
        Some(config) => (Some(config.name), Some(config.id), Some(config.completed)), 
        None => (None, None, None)
    };

    ModelConfig::new(
        sample_type,
        sample_group,
        ercc_input_mass,
        library_tag,
        run_id,
        run_date,
        workflow_name,
        workflow_id, 
        workflow_date
    )
}

#[derive(Deserialize)]
pub struct CerebroInsertSampleQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId
}


#[post("/cerebro")]
async fn insert_model_handler(data: web::Data<AppState>, cerebro: web::Json<Cerebro>, query: web::Query<CerebroInsertSampleQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let sample_name = cerebro.sample.id.clone();
    let sample_tags = cerebro.sample.tags.clone();
    let sample_workflow = cerebro.workflow.id.clone();

    // Check if a model with the same sample sample name, tags and workflow already exists - this
    // allows for same name+tag combination on different workflow and same name and different tag
    // combination on the same workflow (tags must be same order as well)
    match project_collection
        .find_one(doc! { "sample.id": &sample_name, "workflow.id": &sample_workflow, "sample.tags": &sample_tags}, None)
        .await
    {
        Ok(None) => {},
        Ok(Some(_)) => {
            return HttpResponse::Conflict().json(
                serde_json::json!({"status": "fail", "message": format!("Sample {} from workflow {} with the same tags ({}) already exists in database", &sample_name, &sample_workflow, &sample_tags.join(", ")), "data": []})
            )
        }
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": err.to_string(), "data": []})
        ),
    }

    let result = project_collection.insert_one(cerebro.into_inner(), None).await;
    match result {
        Ok(_) => HttpResponse::Ok().json(
            serde_json::json!({"status": "success", "message": format!("Sample {} from workflow {} with tags ({}) added to database", &sample_name, &sample_workflow, &sample_tags.join(", ")), "data": []})
        ),
        Err(err) => HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string()), "data": []})
        ),
    }
}

#[derive(Deserialize)]
struct CerebroGetQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Optional parameters
    taxa: Option<bool>
}

#[get("/cerebro/{id}")]
async fn get_cerebro_uuid(data: web::Data<AppState>, id: web::Path<CerebroId>, query: web::Query<CerebroGetQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_uuid_cerebro_pipeline(&id, &query.taxa);
    
    match project_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let samples = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match samples.is_empty() {
                false => HttpResponse::Ok().json(
                    serde_json::json!({"status": "success", "message": "Documents found", "data": samples})
                ),
                true => HttpResponse::Ok().json(
                    serde_json::json!({"status": "fail", "message": "No document found", "data": []})
                )
            }
        },
        Err(err) => HttpResponse::InternalServerError().body(err.to_string()),
    }
}

#[derive(Deserialize)]
struct CerebroAddPriorityTaxonQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Comma-separated CerebroIds
    id: CerebroIds
}

#[post("/cerebro/priority-taxa")]
async fn add_priority_taxon_handler(request: HttpRequest, data: web::Data<AppState>, priority_taxon: web::Json<PriorityTaxonSchema>, query: web::Query<CerebroAddPriorityTaxonQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let ids = query.id.split(",").map(|x| x.trim()).collect::<Vec<&str>>();
    
    // Disable comments on the priority taxon if global security is set
    // this may not be necessary for front-end, but is done for assurance
    // that no comments (even if submitted as manual request to API) are
    // stored in the database
    let priority_taxon = match !data.env.security.components.comments {
        true => {
            let mut pt = priority_taxon.into_inner();
            pt.comment = String::new();
            PriorityTaxon::from_schema(pt)
        },
        false => {
            PriorityTaxon::from_schema(priority_taxon.into_inner())
        }
    };

    let update = doc! { "$push": { "sample.priority": &mongodb::bson::to_bson(&priority_taxon).unwrap() } };  // unwrap call see if need to handle

    match project_collection
        .update_many(
            doc! { "id":  { "$in" : &ids } }, update, None)
        .await
    {   
        Ok(_) => {
            // Log the action in admin and team databases
            match log_database_change(
                &data, 
                auth_guard.team, 
                RequestLog::new(
                    LogModule::UserAction,
                    Action::PriorityTaxonAdded,
                    false,
                    priority_taxon.log_description(),
                    AccessDetails::new( 
                        &request, 
                        Some(&auth_guard.user.id),
                        Some(&auth_guard.user.email),
                        Some(&query.db), 
                        Some(&query.project)
                    )
                )
            ).await
            {
                Ok(_) => HttpResponse::Ok().json(
                    serde_json::json!({"status": "success", "message": "Added priority taxon for requested documents"})
                ),
                Err(err_response) => err_response
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string())})
        ),
    }
}

// Keep query for DB/Project authorized reports from team collection

#[allow(dead_code)]
#[derive(Deserialize)]
struct CerebroGetStoredReportQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
}

#[get("/cerebro/reports/{id}")]
async fn get_report_from_storage_handler(data: web::Data<AppState>, id: web::Path<String>, _: web::Query<CerebroGetStoredReportQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let reports_collection: Collection<ReportEntry> = get_teams_db_collection(&data, auth_guard.team, TeamAdminCollection::Reports);
    match reports_collection
        .find_one(
            doc! { "id":  &id.into_inner() }, None)
        .await
    {   
        Ok(None) => {
            HttpResponse::Ok().json(
                serde_json::json!({
                    "status": "fail", 
                    "message": "Could not find stored report",
                    "data": serde_json::json!({}) 
                })
            )
        },
        Ok(Some(report_entry)) => {
            HttpResponse::Ok().json(
                serde_json::json!({
                    "status": "success", 
                    "message": "Returned stored report",
                    "data": serde_json::json!({"report_entry": report_entry}) 
                })
            )
        },
        Err(_) => HttpResponse::InternalServerError().json(
            serde_json::json!({
                "status": "error", 
                "message": "Failed to query stored report",
                "data": serde_json::json!({}) 
            })
        ),
    }

}

#[derive(Deserialize)]
struct CerebroDeleteStoredReportQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId
}

#[delete("/cerebro/reports/{id}")]
async fn delete_stored_report_handler(request: HttpRequest, data: web::Data<AppState>, id: web::Path<String>, query: web::Query<CerebroDeleteStoredReportQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to delete a report", "data": serde_json::json!({})
        }))
    }

    // TODO: Delete from DB models

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };
    
    let report_id: String = id.into_inner();

    let ids = match project_collection
        .find(
            doc! { "sample.reports.id":  &report_id }, None)
        .await
    {
        Ok(cursor) => cursor
            .try_collect()
            .await
            .unwrap_or_else(|_| vec![])
            .into_iter()
            .map(|c| c.id )
            .collect::<Vec<String>>(),
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error", "message": format!("{}", err.to_string())})
        ),
    };

    if ids.is_empty() {
        return HttpResponse::NotFound().json(
            serde_json::json!({
                "status": "fail", 
                "message": format!("Failed to find documents with report entries for report identifier: {}", &report_id), 
                "data": serde_json::json!({})
        }))
    };

    let update =  doc! { "$pull": { "sample.reports": {"id": &report_id } } };
    match project_collection
        .update_many(
            doc! { "id":  { "$in" : &ids } }, update, None)
        .await
    {   
        Ok(_) => {},
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string()), "data": serde_json::json!({})})
        ),
    };

    // Delete stored report in database but also the entry in the Cerebro models under the shared identifier
    let reports_collection: Collection<ReportEntry> = get_teams_db_collection(&data, auth_guard.team.clone(), TeamAdminCollection::Reports);

    match reports_collection
        .find_one_and_delete(
            doc! { "id":  &report_id}, None)
        .await
    {   
        Ok(deleted) => {
            match deleted {
                Some(_) | None => {
                    // Log the action in admin and team databases
                    match log_database_change(
                        &data, 
                        auth_guard.team, 
                        RequestLog::new(
                            LogModule::UserAction,
                            Action::ReportEntryRemoved,
                            true,  // should not happen except in dev for archival / data retention policy compliance (deletion can be deactivated by admin)
                            format!("Report entry removed: id={}", &report_id),
                            AccessDetails::new( 
                                &request, 
                                Some(&auth_guard.user.id),
                                Some(&auth_guard.user.email),
                                Some(&query.db), 
                                Some(&query.project)
                            )
                        )
                    ).await
                    {
                        Ok(_) => HttpResponse::Ok().json(
                            serde_json::json!({
                                "status": "success", 
                                "message": "Deleted report",
                                "data": serde_json::json!({}) 
                            })
                        ),
                        Err(_) =>  HttpResponse::InternalServerError().json(
                            serde_json::json!({
                                "status": "fail", 
                                "message": "Failed to log deletion of stored report",
                                "data": serde_json::json!({}) 
                            })
                        )
                    }
                    
                }
            }
            
        },
        Err(_) => HttpResponse::InternalServerError().json(
            serde_json::json!({
                "status": "error", 
                "message": "Failed to query stored report",
                "data": serde_json::json!({}) 
            })
        ),
    }

}

#[derive(Deserialize)]
struct CerebroAddReportPdfQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId
}

#[post("/cerebro/reports/pdf")]
async fn create_pdf_report_handler(request: HttpRequest, data: web::Data<AppState>, report_schema: web::Json<ReportSchema>, query: web::Query<CerebroAddReportPdfQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {


    if !auth_guard.user.roles.contains(&Role::Report) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have sufficient permissions to create a report",
        }))
    }

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let report = match build_report(&project_collection, &report_schema).await {
        Ok(report) => report,
        Err(err_response) => return err_response
    };

    let report_id = report.id.clone();

    let rendered = actix_web::rt::task::spawn_blocking(move || {
        report.render_pdf(None, None, None)
    });

    let (report_text, is_pdf) = match rendered.await {
        Ok(result) => {
            match result {
                Ok(pdf_bytes) => (general_purpose::STANDARD.encode(pdf_bytes), true),
                Err(err) => return HttpResponse::InternalServerError().json(
                    serde_json::json!({"status": "fail", "message": format!("Failed to render LaTeX to PDF: {}", err.to_string()), "data": serde_json::json!({}) })
                )
            }
        },
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("Failed to render LaTeX to PDF: {}", err.to_string()), "data": serde_json::json!({}) })
        )
    };
    
    let access_details = AccessDetails::new( 
        &request, 
        Some(&auth_guard.user.id),
        Some(&auth_guard.user.email),
        Some(&query.db), 
        Some(&query.project)
    );
    
    store_and_log_report(&data, auth_guard.team, &project_collection, &report_schema, report_id, report_text, is_pdf, &access_details).await

    
}


#[derive(Deserialize)]
struct CerebroAddReportTexQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId
}

#[post("/cerebro/reports/tex")]
async fn create_typst_report_handler(request: HttpRequest, data: web::Data<AppState>, report_schema: web::Json<ReportSchema>, query: web::Query<CerebroAddReportTexQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    if !auth_guard.user.roles.contains(&Role::Report) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have sufficient permissions to create a report",
        }))
    }

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let report = match build_report(&project_collection, &report_schema).await {
        Ok(report) => report,
        Err(err_response) => return err_response
    };

    let (report_text, is_pdf) = match report.render_template(None) {
        Ok(tex_str) => (tex_str, false),
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("Failed to render report to LaTeX: {}", err.to_string()), "data": serde_json::json!({}) })
        )
    };
    
    let access_details = AccessDetails::new( 
        &request, 
        Some(&auth_guard.user.id),
        Some(&auth_guard.user.email),
        Some(&query.db), 
        Some(&query.project)
    );
    
    store_and_log_report(&data, auth_guard.team, &project_collection, &report_schema, report.id, report_text, is_pdf, &access_details).await
}


#[derive(Deserialize)]
struct CerebroDeletePriorityTaxonQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Comma-separated CerebroIds
    id: CerebroIds
}

#[delete("/cerebro/priority-taxa/{id}")]
async fn delete_priority_taxon_handler(request: HttpRequest, data: web::Data<AppState>, id: web::Path<String>, query: web::Query<CerebroDeletePriorityTaxonQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };
    
    let priority_taxon_id = &id.into_inner();
    let ids = query.id.split(",").map(|x| x.trim()).collect::<Vec<&str>>();
    
    let update =  doc! { "$pull": { "sample.priority": {"id": &priority_taxon_id } } };

    match project_collection
        .update_many(
            doc! { "id":  { "$in" : &ids } }, update, None)
        .await
    {   
        Ok(_) => {
            // Log the action in admin and team databases
            match log_database_change(
                &data, 
                auth_guard.team, 
                RequestLog::new(
                    LogModule::UserAction,
                    Action::PriorityTaxonRemoved,
                    false,
                    format!("Priority taxon removed: id={}", &priority_taxon_id),
                    AccessDetails::new( 
                        &request, 
                        Some(&auth_guard.user.id),
                        Some(&auth_guard.user.email),
                        Some(&query.db), 
                        Some(&query.project)
                    )
                )
            ).await
            {
                Ok(_) => HttpResponse::Ok().json(
                    serde_json::json!({"status": "success", "message": "Removed priority taxon for requested documents"})
                ),
                Err(err_response) => err_response
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string())})
        ),
    }
}


#[derive(Deserialize)]
struct CerebroPatchPriorityTaxonDecisionQuery {
    // Required for access authorization 
    // in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Comma-separated CerebroIDs on which to perform update
    id: CerebroIds
}

#[patch("/cerebro/priority-taxa/decision")]
async fn modify_priority_taxa_decision_handler(request: HttpRequest, data: web::Data<AppState>, decision_schema: web::Json<PriorityTaxonDecisionSchema>, query: web::Query<CerebroPatchPriorityTaxonDecisionQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    // Disable comments on the priority taxon decisionif global security is set
    // this may not be necessary for front-end, but is done for assurance
    // that no comments (even if submitted as manual request to API) are
    // stored in the database
    let decision_schema = match !data.env.security.components.comments {
        true => {
            let mut ds = decision_schema.into_inner();
            ds.decision_comment = String::new();
            ds
        },
        false => {
            decision_schema.into_inner()
        }
    };

    let ids = query.id.split(",").map(|x| x.trim()).collect::<Vec<&str>>();
    
    // Find matching models/priority-taxon to check if it still exists
    match project_collection
        .find(
            doc! { 
                "id":  { "$in" : &ids }, 
                "sample.priority.id": &decision_schema.id
            },
            None 
        )
        .await
    {   
        Ok(cursor) => {
            let priority_taxa = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match priority_taxa.is_empty() {
                true => return HttpResponse::NotFound().json(
                    serde_json::json!({"status": "fail", "message": "Taxon may have been deleted by the user who submitted it, please refresh the page."})
                ),
                false => {} // continue
            }
        },
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string())})
        ),
    }


    // Array filter of which priority taxon to update
    let taxon_update_remove_options = mongodb::options::UpdateOptions::builder().array_filters(
        vec![doc! { "pt.id": &decision_schema.id }]
    ).build();

    // Update #1: We always remove the user's previous decision on the priority taxon first 
    let remove_decision_update = doc! {
        "$pull": { "sample.priority.$[pt].decisions": { "user_id": &auth_guard.user.id }}
    };
    match project_collection.update_many(
        doc! { "id":  { "$in" : &ids } }, remove_decision_update, taxon_update_remove_options)
    .await
    {   
        Ok(_) => {},
        Err(err) => {
            return HttpResponse::InternalServerError().json(
                serde_json::json!({"status": "fail", "message": format!("{}", err.to_string())})
            )
        }
    }

    let taxon_update_options = mongodb::options::UpdateOptions::builder().array_filters(
        vec![doc! { "pt.id": &decision_schema.id }]
    ).build();

    // Update #2: We then add the updated decision and log it
    let decision = PriorityTaxonDecision::from(&decision_schema, &auth_guard.user);
    let update = doc! { "$push": { "sample.priority.$[pt].decisions": &mongodb::bson::to_bson(&decision).unwrap() } };

    match project_collection
        .update_many(
            doc! { "id":  { "$in" : &ids } }, update, taxon_update_options)
        .await
    {   
        Ok(_) => {
            // Log the action in admin and team databases
            match log_database_change(
                &data, 
                auth_guard.team, 
                RequestLog::new(
                    LogModule::UserDecision,
                    match decision.decision { DecisionType::Accept => Action::PriorityTaxonAccepted, DecisionType::Reject => Action::PriorityTaxonRejected },
                    false,
                    decision.log_description(&decision_schema),
                    AccessDetails::new( 
                        &request, 
                        Some(&auth_guard.user.id),
                        Some(&auth_guard.user.email),
                        Some(&query.db), 
                        Some(&query.project)
                    )
                )
            ).await
            {
                Ok(_) => HttpResponse::Ok().json(
                    serde_json::json!({"status": "success", "message": "Updated priority taxa for requested documents", "data": serde_json::json!({"decision": &decision})})
                ),
                Err(err_response) => err_response
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string())})
        ),
    }
}



#[derive(Deserialize)]
struct TaxaQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Required model identifiers and overview switch
    id: Option<CerebroIds>,
    overview: Option<bool>,
    taxid: Option<String>,
    run_id: Option<String>,
    date_range: Option<String>
}


#[post("/cerebro/taxa")]
async fn filtered_taxa_handler(data: web::Data<AppState>, filter_config: web::Json<TaxonFilterConfig>, query: web::Query<TaxaQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    // Extract and parse optional query parameters
    let ids: Option<Vec<CerebroId>> = query
    .id
    .as_ref()
    .map(|i| i.split(',').map(|x| x.trim().to_string()).collect());

    let run_ids: Option<Vec<String>> = query
    .run_id
    .as_ref()
    .map(|r| r.split(',').map(|x| x.trim().to_string()).collect());

    let date_range: Option<Vec<String>> = query
    .date_range
    .as_ref()
    .map(|d| d.split(',').map(|x| x.trim().to_string()).collect());

    let filter_config = filter_config.into_inner();

    // Build MongoDB pipeline
    let pipeline = get_matched_id_taxa_cerebro_pipeline(ids, date_range, run_ids, None);
    
    // Check if we can move the taxon agregation into the aggregation pipelines, maybe even filters - for now ok

    match project_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let cerebro_taxa = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match cerebro_taxa.is_empty() {
                false => {
                    
                    // `cerebro_taxa` is a Vec<Hashmap<TaxId, Taxon>> for each requested 
                    //  Cerebro document but we need to transform from the BSON format
                    let taxon_maps: Vec<HashMap<String, Taxon>> = cerebro_taxa.iter().map(
                        |taxa| mongodb::bson::from_bson(taxa.into()).map_err(|_| {
                            HttpResponse::InternalServerError().json(serde_json::json!({
                                "status": "error", "message": "Failed to transform retrieved taxa"
                            }))
                        }).unwrap()
                    ).collect();

                    // Aggregation loop of taxa evidence for requested - this may need to be
                    // implemented with a task queue, especially if we do not aggregate
                    // at species level later on (using a taxonomy)
                    let mut aggregated_taxa = HashMap::new();
                    for tax_map in taxon_maps {
                        aggregated_taxa = aggregate(&mut aggregated_taxa, &tax_map)
                    }   

                    // let taxa: Vec<Taxon>= aggregated_taxa.values().cloned().collect();
                    

                    // Applying the rank/evidence filters from the post body - i.e. selected by user in filter interface
                    let taxa: Vec<Taxon> = apply_filters(aggregated_taxa.into_values().collect(), &filter_config);

                    if query.overview.is_some_and(|x| x) {
                        // Summarize the evidence fields into a taxon overview, which is the first layer of the front-end `TaxonomyTable`
                        let taxa_overview: Vec<TaxonOverview> = taxa.iter().map(|taxon| { TaxonOverview::from(taxon) }).collect();
                        
                        HttpResponse::Ok().json(serde_json::json!({
                            "status": "success", "message": "Retrieved aggregated taxa overview with summary of evidence", "data": serde_json::json!({"taxa": taxa_overview})
                        }))
                    } else {
                        if let Some(taxid) = &query.taxid {
                            let taxids = taxid.split(",").collect::<Vec<&str>>();

                            let taxa: Vec<Taxon> = taxa.into_iter().filter(|taxon| taxids.contains(&taxon.taxid.as_str())).collect();

                            HttpResponse::Ok().json(serde_json::json!({
                                "status": "success", "message": "Retrieved aggregated taxon by taxid with evidence", "data": serde_json::json!({"taxa": taxa})
                            }))
                        } else {
                            // Return all taxids
                            HttpResponse::Ok().json(serde_json::json!({
                                "status": "success", "message": "Retrieved aggregated taxa with evidence", "data": serde_json::json!({"taxa": taxa})
                            }))
                        }
                    }
                },
                true => HttpResponse::Ok().json(serde_json::json!({
                    "status": "fail", "message": "Failed to retrieve aggregated taxa with evidence", "data": serde_json::json!({"taxa": []})
                }))
            }
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error in retrieving aggregated taxa with evidence", "data": serde_json::json!({"taxa": []})
        }))
    }

}

#[derive(Deserialize)]
struct TaxaSummaryQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // CSV formatted string of data
    csv: bool
}

use futures::{future::ok, stream::once};

#[post("/cerebro/taxa/summary")]
async fn filtered_taxa_summary_handler(data: web::Data<AppState>, schema: web::Json<TaxaSummarySchema>, query: web::Query<TaxaSummaryQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_taxa_summary_pipeline(&schema.sample_ids, &schema.run_ids, &schema.workflow_ids, &schema.workflow_names);
    
    // Check if we can move the taxon agregation into the aggregation pipelines, maybe even filters - for now ok

    match project_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let cerebro_taxa_summary = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match cerebro_taxa_summary.is_empty() {
                false => {

                    let taxon_data: Vec<TaxaSummaryMongoPipeline> = cerebro_taxa_summary.iter().map(
                        |taxa_summary| mongodb::bson::from_bson(taxa_summary.into()).map_err(|_| {
                            HttpResponse::InternalServerError().json(serde_json::json!({
                                "status": "error", "message": "Failed to transform retrieved taxa summary"
                            }))
                        }).unwrap()
                    ).collect();

                    let taxa_overview_data: Vec<Vec<TaxonSummaryOverview>> = taxon_data.into_iter().map(|taxa_summary| {
                        let taxa: Vec<Taxon> = apply_filters(taxa_summary.taxa.clone().into_values().collect(), &schema.filter_config);
                        let taxa_overview: Vec<TaxonOverview> = taxa.iter().map(|taxon| { TaxonOverview::from(taxon) }).collect();
                        let taxa_overview_summary: Vec<TaxonSummaryOverview> = taxa_overview.iter().map(|taxon_overview|{
                            TaxonSummaryOverview::from_taxon_overview(&taxa_summary, &taxon_overview)
                        }).collect();
                        taxa_overview_summary
                    }).collect();

                    let taxa_summary_overview: Vec<TaxonSummaryOverview> = taxa_overview_data.into_iter().flatten().collect();
                                        
                    // Streaming responses due to size
                    
                    match query.csv {
                        true => {
                            let csv_string = as_csv_string(taxa_summary_overview).unwrap();

                            let response = serde_json::to_string(&serde_json::json!({
                                "status": "succes", "message": "Retrieved taxa summaries for requested samples", "data": serde_json::json!({"taxa_summary": [], "csv": csv_string})
                            })).unwrap().try_into_bytes().unwrap();

                            let body = once(ok::<_, actix_web::Error>(response));

                            HttpResponse::Ok().content_type("application/json").streaming(body)
                        },
                        false => {

                            let response = serde_json::to_string(&serde_json::json!({
                                "status": "succes", "message": "Retrieved taxa summaries for requested samples", "data": serde_json::json!({"taxa_summary": taxa_summary_overview, "csv": ""})
                            })).unwrap().try_into_bytes().unwrap();

                            let body = once(ok::<_, actix_web::Error>(response));

                            HttpResponse::Ok().content_type("application/json").streaming(body)
                        }
                    }

                },
                true => HttpResponse::Ok().json(serde_json::json!({
                    "status": "fail", "message": "Failed to retrieve taxa summaries", "data": serde_json::json!({"taxa_summary": [], "csv": ""})
                }))
            }
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error in retrieving taxa summaries", "data": serde_json::json!({"taxa_summary": [], "csv": ""})
        }))
    }

}


#[derive(Deserialize)]
struct CerebroSamplesOverviewPaginatedQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Required for pagination
    page: i64,
    limit: i64,
    // Optional request parameters
    id: Option<SampleId>,
    notag: Option<String>,
    group: Option<String>,
    run: Option<String>,
    workflow: Option<String>,

}

#[get("/cerebro/samples/overview")]
async fn samples_overview_handler(data: web::Data<AppState>, query: web::Query<CerebroSamplesOverviewPaginatedQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let exclude_tags = match &query.notag {
        Some(value) => value.split(",").map(|x| x.trim()).collect::<Vec<&str>>(),
        None => Vec::new()
    };

    let collation_setting = mongodb::options::Collation::builder().locale("en_US").numeric_ordering(true).build();
    let aggregate_options = mongodb::options::AggregateOptions::builder().collation(collation_setting).build();

    let pipeline = get_paginated_sample_overview_pipeline(&query.page, &query.limit, exclude_tags, &query.id, &query.run, &query.workflow, &query.group);

    match project_collection
        .aggregate(pipeline, aggregate_options)
        .await
    {
        Ok(cursor) => {
            let sample_overview = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            
            match sample_overview.is_empty() {
                false => HttpResponse::Ok().json(serde_json::json!({
                    "status": "success", "message": "Retrieved sample data", "data": serde_json::json!({"sample_overview": sample_overview})
                })),
                true => HttpResponse::Ok().json(serde_json::json!({
                    "status": "fail", "message": "No samples available - has data been uploaded?", "data": serde_json::json!({"sample_overview": []})
                })),
            }
            
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error retrieving sample data", "data": serde_json::json!({"sample_overview": []})
        }))
    }

}

#[derive(Deserialize)]
struct CerebroSampleOverviewIdQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Aggregation pipelines to get specific overviews
    workflow: Option<bool>
}

#[get("/cerebro/samples/overview/{id}")]
async fn sample_overview_id_handler(data: web::Data<AppState>, id: web::Path<String>, query: web::Query<CerebroSampleOverviewIdQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_sample_overview_pipeline(&id, &query.workflow);

    match project_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let sample_overview = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            HttpResponse::Ok().json(serde_json::json!({
                "status": "success", "message": "Retrieved sample overview data", "data": serde_json::json!({"sample_overview": sample_overview})
            }))
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error retrieving sample overview data", "data": serde_json::json!({"sample_overview": []})
        })),
    }
}

#[derive(Deserialize)]
struct CerebroSampleIdQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Aggregation pipelines to get specific overviews
    workflow: Option<WorkflowId>
}

#[get("/cerebro/samples/{id}")]
async fn sample_id_handler(data: web::Data<AppState>, id: web::Path<String>, query: web::Query<CerebroSampleIdQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_sample_cerebro_notaxa_pipeline(&id, &query.workflow);

    match project_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let cerebro = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match cerebro.is_empty() {
                false => HttpResponse::Ok().json(serde_json::json!({
                    "status": "success", "message": "Retrieved models by sample", "data": serde_json::json!({"cerebro": cerebro})
                })),
                true => HttpResponse::Ok().json(serde_json::json!({
                    "status": "fail", "message": "No model found matching requested sample", "data": serde_json::json!({"cerebro": []})
                })),
            }
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error retrieving data models for requested sample", "data": serde_json::json!({"cerebro": []})
        })),
    }
}


#[derive(Deserialize)]
struct CerebroSampleSummaryQcQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Aggregation pipelines to get specific overviews
    workflow: Option<WorkflowId>,
    // Return data as comma-separated string
    csv: Option<bool>,
    // Include ERCC biomass estimates if ERCC input mass is provided (picogram)
    ercc: Option<f64>
}


#[derive(Deserialize)]
struct SummaryAggregateResult {
    id: String,
    quality: QualityControl,
    run: RunConfig,
    sample: SampleConfig,
    workflow: WorkflowConfig
}

#[post("/cerebro/samples/summary/qc")]
async fn sample_qc_summary_handler(data: web::Data<AppState>, samples: web::Json<SampleSummaryQcSchema>, query: web::Query<CerebroSampleSummaryQcQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_samples_cerebro_notaxa_pipeline(&samples.sample_ids, &samples.cerebro_ids, &query.workflow);

    if pipeline.is_empty(){
        return HttpResponse::BadRequest().json(serde_json::json!({
            "status": "fail", "message": "Must specify one of sample or model identifiers", "data": serde_json::json!({"summary": [], "csv": ""})
        }))
    };

    match project_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let matched_cerebro = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match matched_cerebro.is_empty() {
                false => {

                    let mut matched_cerebro_data: Vec<SummaryAggregateResult> = matched_cerebro.iter().map(|doc| mongodb::bson::from_bson(doc.into()).map_err(|_| {
                        HttpResponse::InternalServerError().json(serde_json::json!({
                            "status": "error", "message": "Failed to transform retrieved samples to retrieve quality control summary"
                        }))
                    }).unwrap()).collect();

                    let summaries: Vec<ReadQualityControl> = match matched_cerebro_data.iter_mut().map(|result| -> Result<ReadQualityControl, HttpResponse> {

                        match result.quality.reads.with_model(
                                &result.id, 
                                qc_config_from_model(
                                    Some(result.sample.clone()), 
                                    Some(result.run.clone()),
                                    Some(result.workflow.clone()), 
                                )
                            ){
                                Ok(summary) => Ok(summary),
                                Err(_) => Err(HttpResponse::InternalServerError().json(serde_json::json!({
                                    "status": "error", "message": "Failed to transform retrieved sample data into quality control summaries"
                                })))
                        }
                    }).collect() {
                        Ok(summaries) => summaries,
                        Err(err_response) => return err_response
                    };


                    match query.csv {
                        Some(csv) => {
                            if csv {
                                let csv_string = as_csv_string(summaries).unwrap();

                                HttpResponse::Ok().json(serde_json::json!({
                                    "status": "succes", "message": "Retrieved sample summaries for requested samples", "data": serde_json::json!({"summary": [], "csv": csv_string})
                                }))
                            } else {
                                HttpResponse::Ok().json(serde_json::json!({
                                    "status": "succes", "message": "Retrieved sample summaries for requested samples", "data": serde_json::json!({"summary": summaries, "csv": ""})
                                }))
                            }
                        },
                        None => HttpResponse::Ok().json(serde_json::json!({
                            "status": "succes", "message": "Retrieved sample summaries for requested samples", "data": serde_json::json!({"summary": summaries, "csv": ""})
                        }))
                    }
                } ,
                true => HttpResponse::Ok().json(serde_json::json!({
                    "status": "fail", "message": "No data found for requested samples", "data": serde_json::json!({"summary": [], "csv": ""})
                })),
            }
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error retrieving sample summaries for requested samples", "data": serde_json::json!({"summary": [], "csv": ""})
        })),
    }
}

#[derive(Deserialize)]
struct CerebroSampleDeleteQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Aggregation pipelines to get workflow specific samples
    workflow: Option<WorkflowId>
}


#[derive(Deserialize)]
struct DeleteAggregateResult {
    sample_id: String,
}

#[delete("/cerebro/samples")]
async fn delete_sample_handler(data: web::Data<AppState>, samples: web::Json<SampleDeleteSchema>, query: web::Query<CerebroSampleDeleteQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to delete a sample", "data": serde_json::json!({"cerebro": []})
        }))
    }

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_sample_ids_cerebro_pipeline(&samples.sample_id, &query.workflow);

    match project_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let matched_sample_ids = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match matched_sample_ids.is_empty() {
                false => {
                    let delete_sample_ids: Vec<DeleteAggregateResult> = matched_sample_ids.iter().map(|doc| mongodb::bson::from_bson(doc.into()).map_err(|_| {
                        HttpResponse::InternalServerError().json(serde_json::json!({
                            "status": "error", "message": "Failed to transform retrieved sample identifiers to delete sample",  "data": serde_json::json!({"cerebro": []})
                        }))
                    }).unwrap()).collect();

                    // Delete the queried data models
                    match project_collection.delete_many(doc! { "sample.id": { "$in": delete_sample_ids.iter().map(|x| x.sample_id.clone()).collect::<Vec<String>>() }}, None).await {
                        Ok(_) => {
                            HttpResponse::Ok().json(serde_json::json!({
                                "status": "success", "message": "Deleted requested samples", "data": serde_json::json!({"cerebro": []})
                            }))
                        },
                        Err(_) => {
                            HttpResponse::InternalServerError().json(serde_json::json!({
                                "status": "fail", "message": "Failed to delete requested samples", "data": serde_json::json!({"cerebro": []})
                            }))
                        }
                    }
                },
                true => HttpResponse::Ok().json(serde_json::json!({
                    "status": "fail", "message": "No data models found to delete", "data": serde_json::json!({"cerebro": []})
                })),
            }
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error deleting samples", "data": serde_json::json!({"cerebro": []})
        })),
    }
}

#[derive(Deserialize)]
struct CerebroSampleDescriptionUpdateQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Sample identifier to update in this project
    id: SampleId
}

#[patch("/cerebro/samples/description")]
async fn sample_description_handler(data: web::Data<AppState>, update: web::Json<SampleDescriptionSchema>, query: web::Query<CerebroSampleDescriptionUpdateQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to delete a sample", "data": serde_json::json!({"cerebro": []})
        }))
    }

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };


    match project_collection.update_many(doc! { "sample.id": &query.id }, doc! { "$set": { "sample.description": &update.description, "sample.type": &update.sample_type, "sample.group": &update.sample_group } }, None).await
    {
        Ok(_) => HttpResponse::Ok().json(serde_json::json!({
            "status": "success", "message": "Updated sample description", "data": serde_json::json!({})
        })),
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error updating sample description", "data": serde_json::json!({})
        })),
    }
}


#[derive(Deserialize)]
struct CerebroSampleCommentAddQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Comma-separated CerebroIDs for which to add comments
    id: CerebroIds
}

#[post("/cerebro/samples/comment")]
async fn sample_comment_handler(request: HttpRequest, data: web::Data<AppState>, comment: web::Json<SampleCommentSchema>, query: web::Query<CerebroSampleCommentAddQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let comment_schema = comment.into_inner();
    let comment = SampleComment::from_schema(&comment_schema);    
    let ids = query.id.split(",").map(|x| x.trim()).collect::<Vec<&str>>();

    match project_collection.update_many(doc! { "id":  { "$in" : &ids } }, doc! { "$push": { "sample.comments": &mongodb::bson::to_bson(&comment).unwrap() } }, None).await
    {
        Ok(_) => {
            // Log the action in admin and team databases
            match log_database_change(
                &data, 
                auth_guard.team, 
                RequestLog::new(
                    LogModule::UserComment,
                    Action::SampleCommentAdded,
                    false,
                    comment.log_description(),
                    AccessDetails::new( 
                        &request, 
                        Some(&auth_guard.user.id),
                        Some(&auth_guard.user.email),
                        Some(&query.db), 
                        Some(&query.project)
                    )
                )
            ).await {
                Ok(_) => HttpResponse::Ok().json(
                    serde_json::json!({"status": "success", "message": "Added sample comment for requested documents", "data": serde_json::json!({"comment": []})})
                ),
                Err(err_response) => err_response
            }
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error adding sample comment", "data": serde_json::json!({})
        })),
    }
}

#[derive(Deserialize)]
struct CerebroSampleCommentDeleteQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Comma-separated CerebroIDs for which to add comments
    id: CerebroIds,
    // Comemnt identifier to delete
    comment_id: CerebroIds
}

#[delete("/cerebro/samples/comment")]
async fn delete_sample_comment_handler(request: HttpRequest, data: web::Data<AppState>, query: web::Query<CerebroSampleCommentDeleteQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let ids = query.id.split(",").map(|x| x.trim()).collect::<Vec<&str>>();

    match project_collection.update_many(doc! { "id":  { "$in" : &ids } }, doc! { "$pull": { "sample.comments": { "comment_id": &query.comment_id }} }, None).await
    {
        Ok(_) => {
            // Log the action in admin and team databases
            match log_database_change(
                &data, 
                auth_guard.team, 
                RequestLog::new(
                    LogModule::UserComment,
                    Action::SampleCommentRemoved,
                    false,
                    format!("Sample comment removed: comment_id={}", &query.comment_id),
                    AccessDetails::new( 
                        &request, 
                        Some(&auth_guard.user.id),
                        Some(&auth_guard.user.email),
                        Some(&query.db), 
                        Some(&query.project)
                    )
                )
            ).await {
                Ok(_) => HttpResponse::Ok().json(
                    serde_json::json!({"status": "success", "message": "Removed sample comment for requested documents", "data": serde_json::json!({"comment": []})})
                ),
                Err(err_response) => err_response
            }
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error removing sample comments", "data": serde_json::json!({})
        })),
    }
}

#[derive(Deserialize)]
struct CerebroSampleQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
}

#[get("/cerebro/samples")]
async fn sample_handler(data: web::Data<AppState>, query: web::Query<CerebroSampleQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    match project_collection
    .distinct("sample", None, None) // All distinct samples in this project
    .await
    {
        Ok(samples) => {
            HttpResponse::Ok().json(serde_json::json!({
                "status": "success", "message": "Retrieved distinct sample configurations", "data": serde_json::json!({"samples": samples})
            }))
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error retrieving distinct sample configurations", "data": serde_json::json!({"samples": []})
        }))
    }
}


#[derive(Deserialize)]
struct CerebroWorkflowQuery {
    // Required for access authorization 
    // in user guard middleware
    db: DatabaseId,
    project: ProjectId
}

#[get("/cerebro/workflows")]
async fn workflow_handler(data: web::Data<AppState>,query: web::Query<CerebroWorkflowQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    match project_collection
    .distinct("workflow", None, None) // All distinct workflows in this project
    .await
    {
        Ok(workflows) => {
            HttpResponse::Ok().json(serde_json::json!({
                "status": "success", "message": "Retrieved distinct workflow configurations", "data": serde_json::json!({"workflows": workflows})
            }))
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error retrieving distinct workflow configurations", "data": serde_json::json!({"samples": []})
        }))
    }
}

#[derive(Deserialize)]
struct CerebroRunQuery {
    // Required for access authorization 
    // in user guard middleware
    db: DatabaseId,
    project: ProjectId
}

#[get("/cerebro/runs")]
async fn run_handler(data: web::Data<AppState>, query: web::Query<CerebroRunQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    match project_collection
    .distinct("run", None, None) // All distinct runs in this project
    .await
    {
        Ok(runs) => {
            HttpResponse::Ok().json(serde_json::json!({
                "status": "success", "message": "Retrieved distinct run configurations", "data": serde_json::json!({"runs": runs})
            }))
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error retrieving distinct run configurations", "data": serde_json::json!({"samples": []})
        })),
    }
}
    

#[derive(Deserialize)]
struct CerebroWorkflowIdQuery {
    // Required for access authorization 
    // in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Aggregation pipeline options
    runs: Option<String>,
    tags: Option<String>
}

#[get("/cerebro/workflows/{id}")]
async fn workflow_id_handler(data: web::Data<AppState>, id: web::Path<WorkflowId>, query: web::Query<CerebroWorkflowIdQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_workflow_cerebro_notaxa_pipeline(&id, &query.runs, &query.tags, &None);

    match project_collection
        .aggregate(pipeline, None)
        .await
    {
        Ok(cursor) => {
            let cerebro = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match cerebro.is_empty() {
                false => HttpResponse::Ok().json(serde_json::json!({
                    "status": "success", "message": "Retrieved models by workflow", "data": serde_json::json!({"cerebro": cerebro})
                })),
                true => HttpResponse::NotFound().json(serde_json::json!({
                    "status": "fail", "message": "Failed to retrieve models by workflow", "data": serde_json::json!({"cerebro": []})})
                ),
            }
        },
        Err(_) =>  HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error retrieving models by workflow", "data": serde_json::json!({"cerebro": []})
        })),
    }
}

#[get("/cerebro/status")]
async fn status_handler(_: web::Data<AppState>, _: jwt::JwtDataMiddleware) -> HttpResponse {
    HttpResponse::Ok().into()
}

// ================
// Helper functions
// ================

type MongoDatabaseName = String;

pub fn get_authorized_database_and_project_collection(data: &web::Data<AppState>, db: &DatabaseId, project: &ProjectId, user_guard: &jwt::JwtDataMiddleware) -> Result<(MongoDatabaseName, Collection<Cerebro>), HttpResponse> {
    
    let (database, project) = {
            let database_matches: Vec<&TeamDatabase> = user_guard.team.databases.iter().filter(|x| &x.id == db).collect();

            if database_matches.len() != 1 {
                return Err(HttpResponse::NotFound().json(serde_json::json!({"status": "fail", "message": "Could not find requested database"})))
            }

            let database_match = database_matches[0];
            let project_matches: Vec<&ProjectCollection> = database_match.projects.iter().filter(|x| &x.id == project).collect();

            if project_matches.len() != 1 {
                return Err(HttpResponse::NotFound().json(serde_json::json!({"status": "fail", "message": "Could not find requested project"})))
            }

            let project_match = project_matches[0];

            (database_match.to_owned(), project_match.to_owned())
    };

    let mongo_db = data.db.database(&database.database);
    Ok((database.database, mongo_db.collection(&project.collection)))
}


async fn build_report(project_collection: &Collection<Cerebro>, report_schema: &web::Json<ReportSchema>) -> Result<ClinicalReport, HttpResponse> {

    // Find matching models
    let mut cerebro = match project_collection
        .find(
            doc! { 
                "id":  { "$in" : &report_schema.ids }
            },
            None 
        )
        .await
    {   
        Ok(cursor) => {
            let cerebro = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match cerebro.is_empty() {
                true => return Err(HttpResponse::Ok().json(
                    serde_json::json!({"status": "fail", "message": "Failed to find requested models for reporting" })
                )),
                false => {
                    if cerebro.len() != report_schema.ids.len() {
                        return Err(HttpResponse::Ok().json(
                            serde_json::json!({"status": "fail", "message": "Failed to find the exact requested models for reporting"})
                        ))
                    }
                    cerebro
                }
            }
        },
        Err(err) => return Err(HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string())})
        )),
    };

    // Compute quality control summaries from the models
    let quality_summaries: Vec<ReadQualityControl> = match cerebro.iter_mut().map(|result| -> Result<ReadQualityControl, HttpResponse> {

        match result.quality.reads.with_model(
            &result.id,
            qc_config_from_model(
                Some(result.sample.clone()), 
                Some(result.run.clone()),
                Some(result.workflow.clone()), 
            )
         ) {
            Ok(summary) => Ok(summary),
            Err(_) => Err(HttpResponse::InternalServerError().json(serde_json::json!({
                "status": "error", "message": "Failed to transform retrieved sample data into quality control summaries"
            })))
        }
    }).collect() {
        Ok(summaries) => summaries,
        Err(err_response) => return Err(err_response)
    };

    // Read the requested template from file
    let template_str = match std::fs::read_to_string("/data/templates/report/template.toml") {
        Ok(template_str) => template_str,
        Err(_) => return Err(HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to read the report template configuration file" })
        ))
    };

    // Read the requested template from str
    let mut template_config: TemplateConfig = match toml::from_str(&template_str) {
        Ok(config) => config,
        Err(_) => return Err(HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to read the report template configuration TOML" })
        ))
    }; 

    // Must be set to find path to the underlying TeX Handlebars template file and PNG logo for render
    template_config.file_path = PathBuf::from("/data/templates/report/report.hbs");
    template_config.logo_path = PathBuf::from("/data/templates/report/logo.png");

    // Read the requested template from file
    let assay_str = match std::fs::read_to_string("/data/templates/report/assay.toml") {
        Ok(assay_str) => assay_str,
        Err(_) => return Err(HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to read the assay template configuration file" })
        ))
    };

    // Read the requested assay template from file
    let assay_config: AssayTemplate = match toml::from_str(&assay_str) {
        Ok(config) => config,
        Err(_) => return Err(HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to read the assay template configuration TOML" })
        ))
    }; 

     // Read the requested template from file
     let bioinfo_str = match std::fs::read_to_string("/data/templates/report/bioinformatics.toml") {
        Ok(bioinfo_str) => bioinfo_str,
        Err(_) => return Err(HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to read the bioinformatics template configuration file" })
        ))
    };

    // Read the requested assay template from file
    let bioinfo_config: BioinformaticsTemplate = match toml::from_str(&bioinfo_str) {
        Ok(config) => config,
        Err(_) => return Err(HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to read the bioinformatics template configuration TOML" })
        ))
    }; 

    // Create the report from schema and templates
    Ok(
        ClinicalReport::from_api(&report_schema, &mut template_config, &assay_config, &bioinfo_config, &quality_summaries)
    )
    

}

async fn store_and_log_report(
    data: &web::Data<AppState>,
    team: Team,
    project_collection: &Collection<Cerebro>,  
    report_schema: &web::Json<ReportSchema>,
    report_id: uuid::Uuid,
    report_text: String,
    is_pdf: bool,
    access_details: &AccessDetails
) -> HttpResponse {

    // Create a report entry for team database storage that includes the compiled report
    let reports_collection: Collection<ReportEntry> = get_teams_db_collection(&data, team.clone(), TeamAdminCollection::Reports);
    let mut report_entry = ReportEntry::from_schema(report_id.to_string(), &report_schema, Some(report_text.clone()), Some(is_pdf));

    match reports_collection.insert_one(
        &report_entry, None
    ).await {
        Ok(_) => {},
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to insert report into storage" })
        )
    }

    // Modify report entry to maintain same UUID and to log the report creation in the requested models themselves
    // this is done without the text to avoid overloading model sizes when many reports are generated
    report_entry.report_text = None;
    report_entry.report_pdf = None;

    let update = doc! { "$push": { "sample.reports": &mongodb::bson::to_bson(&report_entry).unwrap() } };  // unwrap call see if need to handle

    match project_collection
        .update_many(
            doc! { "id":  { "$in" : &report_schema.ids } }, update, None)
        .await
    {   
        Ok(_) => {
            // Log the action in admin and team databases
            match log_database_change(
                &data, 
                team, 
                RequestLog::new(
                    LogModule::UserAction,
                    Action::ReportEntryAdded,
                    false,
                    report_entry.log_description(),
                    access_details.clone(),
                )
            ).await
            {
                Ok(_) => {
                    HttpResponse::Ok().json(
                        serde_json::json!({
                            "status": "success", 
                            "message": "Generated report (PDF)", 
                            "data": serde_json::json!({"report": report_text, "pdf": is_pdf}) 
                        })
                    )
                },
                Err(err_response) => err_response
            }
        },
        Err(_) => HttpResponse::InternalServerError().json(
            serde_json::json!({
                "status": "fail", 
                "message": "Failed to update models with report entry"
            })
        ),
    }
}

pub fn cerebro_config(cfg: &mut web::ServiceConfig, config: &Config) {
    
    cfg.service(insert_model_handler)
       .service(samples_overview_handler)
       .service(sample_overview_id_handler)
       .service(sample_id_handler)
       .service(workflow_id_handler)
       .service(filtered_taxa_handler)
       .service(add_priority_taxon_handler)
       .service(delete_priority_taxon_handler)
       .service(modify_priority_taxa_decision_handler)
       .service(workflow_handler)
       .service(run_handler)
       .service(sample_handler)
       .service(delete_sample_handler)
       .service(sample_qc_summary_handler)
       .service(sample_description_handler)
       .service(filtered_taxa_summary_handler)
       .service(create_typst_report_handler)
       .service(get_report_from_storage_handler)
       .service(delete_stored_report_handler)
       .service(create_pdf_report_handler)
       .service(status_handler);

    if config.security.components.comments {
        // Disable comments on sample and priority taxon decision comment modification endpoints
        // Comments on priority taxon decisions and submissions are disabled in the endpoints
        cfg.service(sample_comment_handler)
           .service(delete_sample_comment_handler);
    }

    

}