use cerebro_model::api::cerebro::response::{CerebroIdentifierResponse, CerebroIdentifierSummary, ContaminationTaxaResponse, FilteredTaxaResponse, SampleSummaryResponse, SampleSummary, TaxonHistoryResponse, TaxonHistoryResult};
use futures::StreamExt;
use mongodb::bson::{from_document, Document};
use mongodb::gridfs::GridFsBucket;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use futures::stream::TryStreamExt;
use mongodb::{bson::doc, Collection, Database};
use actix_web::{get, post, web, HttpResponse, patch, delete, HttpRequest};

use crate::api::auth::jwt::{self, TeamProjectAccessQuery};
use cerebro_model::api::config::Config;
use cerebro_model::api::users::model::Role;
use crate::api::utils::as_csv_string;
use crate::api::server::AppState;
use crate::api::cerebro::{gridfs, mongo::*};
use crate::api::logs::utils::log_database_change;
use cerebro_model::api::logs::model::{LogModule, Action, RequestLog, AccessDetails};
use cerebro_model::api::teams::model::{DatabaseId, ProjectCollection, ProjectId, Team, TeamAdminCollection, TeamDatabase};
use cerebro_model::api::cerebro::model::{
    Cerebro, 
    PriorityTaxon,
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
    CerebroIdentifierSchema, ContaminationSchema, PriorityTaxonDecisionSchema, PriorityTaxonSchema, SampleCommentSchema, SampleDeleteSchema, SampleDescriptionSchema, SampleGroupSchema, SampleSummarySchema, SampleTypeSchema, TaxaSummarySchema
};

use mongodb::options::{UpdateManyModel, WriteModel};

use cerebro_pipeline::taxa::filter::*;
use cerebro_pipeline::modules::quality::{QualityControl, ReadQualityControl, ModelConfig};
use cerebro_pipeline::taxa::taxon::{aggregate, collapse_taxa, Taxon};

type CerebroIds = String;

fn qc_config_from_model(sample_config: Option<SampleConfig>, run_config: Option<RunConfig>, workflow_config: Option<WorkflowConfig>) -> ModelConfig {

    let (sample_type, sample_group, ercc_input_mass, library_tag) = match sample_config {
        Some(config) => (config.sample_type, config.sample_group, config.ercc_input_mass, Some(config.tags.join("-"))), 
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


#[post("/cerebro")]
async fn insert_model_handler(data: web::Data<AppState>, cerebro: web::Json<Cerebro>, auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let (_, db, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
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
        .find_one(doc! { "sample.id": &sample_name, "workflow.id": &sample_workflow, "sample.tags": &sample_tags})
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

    let mut model = cerebro.into_inner();

    match gridfs::upload_taxa_to_gridfs(db.gridfs_bucket(None), &model.taxa, &model.id).await {
        Ok(_) => model.taxa.clear(), // Set taxa to empty since we uploaded to GridFS
        Err(err) => {
            return HttpResponse::InternalServerError().json(serde_json::json!({
                "status": "fail",
                "message": format!("GridFS upload error: {}", err),
                "data": []
            }));
        }
    }

    let result = project_collection.insert_one(model).await;

    match result {
        Ok(_) => HttpResponse::Ok().json(
            serde_json::json!({"status": "success", "message": format!("Sample {} from workflow {} with tags ({}) added to database", &sample_name, &sample_workflow, &sample_tags.join(", ")), "data": []})
        ),
        Err(err) => HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("{}", err.to_string()), "data": []})
        ),
    }
}

#[get("/cerebro")]
async fn retrieve_model_handler(data: web::Data<AppState>, auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let (_, db, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    match project_collection
        .find(doc! {})
        .await
    {
        Ok(cursor) => {

            let mut samples = cursor
                .try_collect()
                .await
                .unwrap_or_else(|_| vec![]);

            match samples.is_empty() {
                false => {

                    let mut models = Vec::new();
                    for model in samples.iter_mut() {
                        
                        match gridfs::download_taxa_from_gridfs(db.gridfs_bucket(None), &model.id).await {
                            Ok(taxa) => {
                                model.taxa = taxa;
                                models.push(model)
                            }
                            Err(err) => {
                                return HttpResponse::InternalServerError().json(format!("GridFS download error: {}", err));
                            }
                        };
                    }
                    
                    HttpResponse::Ok().json(
                        serde_json::json!({"status": "success", "message": "Documents found", "data": models})
                    )
                },
                true => HttpResponse::Ok().json(
                    serde_json::json!({"status": "fail", "message": "No document found", "data": []})
                )
            }
        },
        Err(err) => HttpResponse::InternalServerError().body(err.to_string()),
    }
}

#[derive(Deserialize)]
struct CerebroGetQuery {
    // Optional parameters
    taxa: Option<bool>
}

#[get("/cerebro/{id}")]
async fn get_cerebro_uuid(data: web::Data<AppState>, id: web::Path<CerebroId>, query: web::Query<CerebroGetQuery>, auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, db, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_uuid_cerebro_pipeline(&id, &query.taxa);
    
    match project_collection
        .aggregate(pipeline)
        .await
    {
        Ok(cursor) => {
            let samples = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match samples.is_empty() {
                false => {

                    let mut models = Vec::new();
                    for doc in samples {

                        let mut model: Cerebro = match mongodb::bson::from_document(doc) {
                            Ok(m) => m,
                            Err(err) => {
                                return HttpResponse::InternalServerError().json(
                                    serde_json::json!({
                                        "status": "fail",
                                        "message": format!("Failed to convert document to Cerebro model: {}", err),
                                        "data": []
                                    })
                                );
                            }
                        };
                        
                        match gridfs::download_taxa_from_gridfs(db.gridfs_bucket(None), &model.id).await {
                            Ok(taxa) => {
                                model.taxa = taxa;
                                models.push(model)
                            }
                            Err(err) => {
                                return HttpResponse::InternalServerError().json(format!("GridFS download error: {}", err));
                            }
                        };
                    }
                    
                    HttpResponse::Ok().json(
                        serde_json::json!({"status": "success", "message": "Documents found", "data": models})
                    )
                },
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
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
            doc! { "id":  { "$in" : &ids } }, update)
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };
    
    let priority_taxon_id = &id.into_inner();
    let ids = query.id.split(",").map(|x| x.trim()).collect::<Vec<&str>>();
    
    let update =  doc! { "$pull": { "sample.priority": {"id": &priority_taxon_id } } };

    match project_collection
        .update_many(
            doc! { "id":  { "$in" : &ids } }, update)
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
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
            }
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
    // Build an update model to remove the user's previous decision.
    let remove_model = WriteModel::UpdateMany(
        UpdateManyModel::builder()
            .namespace(project_collection.namespace())
            .filter(doc! { "id": { "$in": &ids } })
            .update(doc! { 
                "$pull": { "sample.priority.$[pt].decisions": { "user_id": &auth_guard.user.id } } 
            })
            .array_filters(Some(vec![mongodb::bson::Bson::Document(doc! { "pt.id": &decision_schema.id })]))
            .build()
    );

    // Build an update model to add the new decision.
    let decision = PriorityTaxonDecision::from(&decision_schema, &auth_guard.user);
    let add_model = WriteModel::UpdateMany(
        UpdateManyModel::builder()
            .namespace(project_collection.namespace())
            .filter(doc! { "id": { "$in": &ids } })
            .update(doc! { 
                "$push": { "sample.priority.$[pt].decisions": mongodb::bson::to_bson(&decision).unwrap() } 
            })
            .array_filters(Some(vec![mongodb::bson::Bson::Document(doc! { "pt.id": &decision_schema.id })]))
            .build()
    );

    // Execute both update operations in a single bulk write.
    match data.db.bulk_write(vec![remove_model, add_model]).await
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


#[post("/cerebro/ids")]
async fn cerebro_identifier_handler(data: web::Data<AppState>, body: web::Json<CerebroIdentifierSchema>, auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    // Primary sample request check
    match project_collection
        .find_one(doc! { "sample.id": &body.sample })
        .await {
            Ok(Some(_)) => {},
            Ok(None) => return HttpResponse::NotFound().json(
                CerebroIdentifierResponse::sample_not_found(&body.sample)
            ),
            Err(err) => return HttpResponse::InternalServerError().json(
                CerebroIdentifierResponse::server_error(err.to_string())
            )
        }

    // Build the match conditions dynamically.
    let mut match_conditions = vec![];

    // If controls are provided, match primary sample or controls; otherwise, just the primary sample.
    if let Some(controls) = &body.controls {
        match_conditions.push(doc! {
            "$or": [
                { "sample.id": &body.sample },
                { "sample.id": { "$in": controls } }
            ]
        });
    } else {
        match_conditions.push(doc! {
            "sample.id": &body.sample
        });
    }

    // If tags are provided, add a condition to filter on sample tags.
    if let Some(tags) = &body.tags {
        match_conditions.push(doc! {
            "sample.tags": { "$in": tags }
        });
    }

    let pipeline = vec![
        // Combine the dynamically built match conditions.
        doc! { "$match": { "$and": match_conditions } },
        // Project only the fields required for CerebroIdentifierSummary.
        doc! { "$project": {
            "_id": 0,
            "cerebro_id": "$id",
            "sample_id": "$sample.id",
            "sample_tags": "$sample.tags",
            "sample_type": "$sample.sample_type",
            "run_id": "$run.id",
            "run_date": "$run.date",
            "workflow_id": "$workflow.id",
            "workflow_name": "$workflow.name"
        } },
    ];

    match project_collection
        .aggregate(pipeline)
        .await
        {
            Ok(cursor) => {
                
                let watchers: Vec<CerebroIdentifierSummary> = cursor
                    .try_collect::<Vec<_>>()
                    .await
                    .unwrap_or_else(|_| vec![])
                    .into_iter()
                    .filter_map(|doc| from_document(doc).ok())
                    .collect();
    
                match watchers.is_empty() {
                    false => HttpResponse::Ok().json(
                        CerebroIdentifierResponse::success(watchers)
                    ),
                    true => HttpResponse::NotFound().json(
                        CerebroIdentifierResponse::not_found()
                    )
                }
            },
            Err(err) => HttpResponse::InternalServerError().json(
                CerebroIdentifierResponse::server_error(err.to_string())
            )
        }
}

#[derive(Deserialize)]
struct TaxaQuery {
    // Required model identifiers and overview switch
    id: Option<CerebroIds>,
    overview: Option<bool>,
    taxid: Option<String>,
    run_id: Option<String>,
    date_range: Option<String>,
    no_evidence: Option<bool>
}


#[derive(Deserialize, Debug)]
struct TaggedTaxa {
    pub id: String,
    pub taxa: Vec<Taxon>,
    pub name: String,
    pub sample_tags: Vec<String>
}

#[post("/cerebro/taxa")]
async fn filtered_taxa_handler(data: web::Data<AppState>, filter_config: web::Json<TaxonFilterConfig>, query: web::Query<TaxaQuery>, auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let (_, db, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
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

    let pipeline = get_matched_id_taxa_cerebro_pipeline(ids, date_range, run_ids);
    
    // Check if we can move the taxon agregation into the aggregation pipelines, maybe even filters - for now ok

    match project_collection
        .aggregate(pipeline)
        .await
    {
        Ok(cursor) => {
            let docs = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
            match docs.is_empty() {
                false => {
                    
                    let mut tagged_taxa = Vec::new();
                    for doc in docs {

                        let mut model: TaggedTaxa = match mongodb::bson::from_document(doc) {
                            Ok(m) => m,
                            Err(err) => {
                                return HttpResponse::InternalServerError().json(
                                    serde_json::json!({
                                        "status": "fail",
                                        "message": format!("Failed to convert document to TaggedTaxa model: {}", err),
                                        "data": []
                                    })
                                );
                            }
                        };
                        
                        match gridfs::download_taxa_from_gridfs(db.gridfs_bucket(None), &model.id).await {
                            Ok(taxa) => {
                                model.taxa = taxa;
                                tagged_taxa.push(model)
                            }
                            Err(err) => {
                                return HttpResponse::InternalServerError().json(format!("GridFS download error: {}", err));
                            }
                        };
                    }

                    // Aggregation loop of taxa evidence for requested samples
                    let mut aggregated_taxa = HashMap::new();
                    for tt in &tagged_taxa {
                        aggregated_taxa = aggregate(&mut aggregated_taxa, &tt.taxa)
                    }               

                    // Create sample name and tags map for filter
                    let mut tag_map = HashMap::new();       
                    for tt in tagged_taxa {
                        tag_map.insert(tt.name, tt.sample_tags);
                    }          
                         

                    // Applying the rank/evidence filters from the provided filter configuration
                    let taxa: Vec<Taxon> = apply_filters(
                        aggregated_taxa.into_values().collect(), 
                        &filter_config, 
                        &tag_map, 
                        match query.no_evidence { Some(no_evidence) => no_evidence, None => false }
                    );
                    
                    if query.overview.is_some_and(|x| x) {                        
                        HttpResponse::Ok().json(FilteredTaxaResponse::success(taxa))
                    } else {
                        // Additional filter for specific taxa by 'taxid':

                        if let Some(taxid) = &query.taxid {

                            let taxids = taxid.split(",").collect::<Vec<&str>>();
                            let taxa: Vec<Taxon> = taxa.into_iter().filter(|taxon| taxids.contains(&taxon.taxid.as_str())).collect();

                            HttpResponse::Ok().json(FilteredTaxaResponse::success(taxa))
                        } else {
                            // Return all otherwise - same as overview at the moment!
                            HttpResponse::Ok().json(FilteredTaxaResponse::success(taxa))
                        }
                    }
                },
                true => HttpResponse::NotFound().json(FilteredTaxaResponse::not_found())
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(FilteredTaxaResponse::server_error(err.to_string()))
    }

}

pub fn apply_prevalence_contamination_filter(
    cerebros: Vec<Cerebro>,
    min_rpm: f64,              // Evidence from any read profiling tool > RPM
    threshold: f64,            // Prevalence percentage threshold (e.g., 0.5 for 50%)
    taxid_subset: Vec<String>, // Only consider these Taxon.taxid values
    collapse_variants: bool
) -> Vec<String> {

    // Check for edge cases
    if cerebros.is_empty() {
        return Vec::new(); // No data to process
    }

    // Convert taxid_subset to a HashSet for efficient lookups
    let taxid_subset_set: HashSet<String> = taxid_subset.into_iter().collect();

    // Total number of Cerebro objects
    let total_cerebro_count = cerebros.len() as f64;

    // HashMap to count how many Cerebro objects have matching Taxon evidence
    let mut taxon_prevalence: HashMap<String, usize> = HashMap::new();

    // Iterate through each Cerebro object
    for cerebro in cerebros {

        // Collect Taxon IDs with matching evidence in this Cerebro object
        let mut found_taxa: HashSet<String> = HashSet::new();

        // If collapse variants is active, collapse the taxon variants (from GTDB)
        // Note this will assign a new taxid to collapsed taxa, which should match
        // any created previously through the 'apply_filter' function

        let cerebro_taxa = if collapse_variants {
           collapse_taxa(cerebro.taxa).expect("FAILED TO COLLAPSE TAXA")
        } else {
            cerebro.taxa
        };

        for taxon in cerebro_taxa {

            // Check if the Taxon is in the subset
            if !taxid_subset_set.is_empty() && !taxid_subset_set.contains(&taxon.taxid) {
                continue; // Skip taxa not in the subset if any are defined in subset, otherwise use Taxon
            }


            // Check if the Taxon has any ProfileRecord evidence matching the criteria
            let has_matching_evidence = taxon
                .evidence
                .profile
                .iter()
                .any(|record| record.rpm >= min_rpm);

            if has_matching_evidence {
                found_taxa.insert(taxon.taxid.clone());
            }
        }

        // Increment prevalence count for each Taxon ID found in this Cerebro
        for taxid in found_taxa {
            *taxon_prevalence.entry(taxid).or_insert(0) += 1;
        }
    }

    // Filter Taxon IDs based on the percentage threshold
    taxon_prevalence
        .into_iter()
        .filter_map(|(taxid, count)| {
            let prevalence = count as f64 / total_cerebro_count;
            if prevalence >= threshold {
                Some(taxid)
            } else {
                None
            }
        })
        .collect()
}

// Prevalence contamination identification
#[post("/cerebro/taxa/contamination")]
async fn contamination_taxa_handler_project(data: web::Data<AppState>, contam_schema: web::Json<ContaminationSchema>, auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let (_, db, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };


    // Obtain all models in this collection with the provided tags
    match project_collection.find(
            doc! { 
                "sample.tags":  { "$in" : &contam_schema.tags }
            }
        ).await {   
        Ok(cursor) => {
            let models = cursor.try_collect().await.unwrap_or_else(|_| vec![]);

            let mut cerebro = Vec::new();
            for mut model in models {
                
                match gridfs::download_taxa_from_gridfs(db.gridfs_bucket(None), &model.id).await {
                    Ok(taxa) => {
                        model.taxa = taxa;
                        cerebro.push(model)
                    }
                    Err(err) => {
                        return HttpResponse::InternalServerError().json(format!("GridFS download error: {}", err));
                    }
                };
            }

            let taxids = apply_prevalence_contamination_filter(cerebro, contam_schema.min_rpm, contam_schema.threshold, contam_schema.taxid.clone(), contam_schema.collapse_variants);
            
            HttpResponse::Ok().json(ContaminationTaxaResponse::success(taxids))
        },
        Err(err) => HttpResponse::InternalServerError().json(ContaminationTaxaResponse::server_error(err.to_string()))
    }
}


#[derive(Deserialize)]
struct TaxaHistoryQuery {
    taxon_label: String,
    host_label: String
}


#[get("/cerebro/taxa/history")]
async fn history_taxa_handler(data: web::Data<AppState>, query: web::Query<TaxaHistoryQuery>, auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {
    
    let (_, db, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    // Match documents that have the provided tax label in their tax_labels array.
    let pipeline = vec![
        doc! { "$match": { "tax_labels": query.taxon_label.clone() } },
        doc! { "$project": {
            "_id": 0,
            "id": "$id",
            "sample_id": "$sample.id",
            "sample_tags": "$sample.tags",
            "sample_type": "$sample.sample_type",
            "run_id": "$run.id",
            "run_date": "$run.date",
            "input_reads": "$quality.reads.input_reads",
            "host_reads": "$quality.reads.host_reads",
            "taxa": "$taxa"
        } },
    ];

    match project_collection.aggregate(pipeline).await {
        Ok(cursor) => {
            
            let docs = cursor.try_collect().await.unwrap_or_else(|_| vec![]);

            match docs.is_empty() {
                false => {
                    
                    let mut history_results: Vec<TaxonHistoryResult> = Vec::new();
                    for doc in docs {
                        let mut model: TaxonHistoryResult = match mongodb::bson::from_document(doc) {
                            Ok(m) => m,
                            Err(err) => {
                                return HttpResponse::InternalServerError().json(
                                    TaxonHistoryResponse::server_error(err.to_string())
                                );
                            }
                        };
                        
                        match gridfs::download_taxa_from_gridfs(db.gridfs_bucket(None), &model.id).await {
                            Ok(taxa) => {
                                model.taxa = taxa.into_iter().filter(|taxon| taxon.lineage.contains(&query.taxon_label) || taxon.lineage.contains(&query.host_label) ).collect();
                                history_results.push(model)
                            }
                            Err(err) => {
                                return HttpResponse::InternalServerError().json(TaxonHistoryResponse::server_error(err.to_string()));
                            }
                        };
                    }

                    HttpResponse::Ok().json(TaxonHistoryResponse::success(history_results))

                },
                true => HttpResponse::NotFound().json(TaxonHistoryResponse::not_found())
            }
        },
        Err(err) => HttpResponse::InternalServerError().json(TaxonHistoryResponse::server_error(err.to_string()))
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let exclude_tags = match &query.notag {
        Some(value) => value.split(",").map(|x| x.trim()).collect::<Vec<&str>>(),
        None => Vec::new()
    };

    let pipeline = get_paginated_sample_overview_pipeline(&query.page, &query.limit, exclude_tags, &query.id, &query.run, &query.workflow, &query.group);

    match project_collection
        .aggregate(pipeline)
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_sample_overview_pipeline(&id, &query.workflow);

    match project_collection
        .aggregate(pipeline)
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
    // Specify workflow to limit the sample identifier query
    workflow: Option<WorkflowId>
}

#[get("/cerebro/samples/{id}")]
async fn sample_id_handler(data: web::Data<AppState>, id: web::Path<String>, query: web::Query<CerebroSampleIdQuery>,  auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_sample_cerebro_notaxa_pipeline(&id, &query.workflow);

    match project_collection
        .aggregate(pipeline)
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


#[derive(Deserialize, Serialize)]
struct CerebroSampleSummaryQuery {
    // Aggregation pipelines to get specific overviews
    workflow: Option<WorkflowId>,
    // Return data in comma-separated table format (string)
    csv: Option<bool>,
    // Include ERCC biomass estimates if ERCC input mass is provided (picogram)
    ercc: Option<f64>
}

#[post("/cerebro/samples/summary/qc")]
async fn sample_qc_summary_handler(data: web::Data<AppState>, samples: web::Json<SampleSummarySchema>, query: web::Query<CerebroSampleSummaryQuery>, auth_query: web::Query<TeamProjectAccessQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &auth_query.db, &auth_query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_samples_cerebro_notaxa_pipeline(&samples.sample_ids, &samples.cerebro_ids, &query.workflow);

    if pipeline.is_empty() {
        return HttpResponse::BadRequest().json(
            SampleSummaryResponse::server_error("Must specify one of sample or model identifiers".to_string())
        );
    }

    match project_collection.aggregate(pipeline).await {
        Ok(cursor) => {
            let matched_docs: Vec<Document> = match cursor.try_collect().await {
                Ok(docs) => docs,
                Err(e) => {
                    return HttpResponse::InternalServerError().json(
                        SampleSummaryResponse::server_error(e.to_string())
                    );
                }
            };

            if matched_docs.is_empty() {
                return HttpResponse::Ok().json(SampleSummaryResponse::not_found());
            }

            let mut matched_data: Vec<SampleSummary> = Vec::new();
            for doc in matched_docs {
                let result: SampleSummary = match mongodb::bson::from_bson(doc.into()) {
                    Ok(res) => res,
                    Err(e) => {
                        return HttpResponse::InternalServerError().json(
                            SampleSummaryResponse::server_error(e.to_string())
                        );
                    }
                };
                matched_data.push(result);
            }

            // Optionally convert the summaries to CSV table (string) if requested.
            let csv_output = if let Some(true) = query.csv {
                match as_csv_string(&matched_data) {
                    Ok(s) => s,
                    Err(response) => return response
                }
            } else {
                "".to_string()
            };

            HttpResponse::Ok().json(SampleSummaryResponse::success(matched_data, csv_output))
        },
        Err(e) => HttpResponse::InternalServerError().json(
            SampleSummaryResponse::server_error(e.to_string())
        ),
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_sample_ids_cerebro_pipeline(&samples.sample_id, &query.workflow);

    match project_collection
        .aggregate(pipeline)
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
                    match project_collection.delete_many(doc! { "sample.id": { "$in": delete_sample_ids.iter().map(|x| x.sample_id.clone()).collect::<Vec<String>>() }}).await {
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
            "status": "fail", "message": "You do not have permission to modify a sample", "data": serde_json::json!({"cerebro": []})
        }))
    }

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };


    match project_collection.update_many(doc! { "sample.id": &query.id }, doc! { "$set": { "sample.description": &update.description } }).await
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
struct CerebroSampleTypeUpdateQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Sample identifier to update in this project
    id: SampleId
}

#[patch("/cerebro/samples/type")]
async fn sample_type_handler(data: web::Data<AppState>, update: web::Json<SampleTypeSchema>, query: web::Query<CerebroSampleTypeUpdateQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to modify a sample", "data": serde_json::json!({"cerebro": []})
        }))
    }

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };


    match project_collection.update_many(doc! { "sample.id": &query.id }, doc! { "$set": { "sample.sample_type": &update.sample_type } }).await
    {
        Ok(_) => HttpResponse::Ok().json(serde_json::json!({
            "status": "success", "message": "Updated sample type", "data": serde_json::json!({})
        })),
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error updating sample type", "data": serde_json::json!({})
        })),
    }
}


#[derive(Deserialize)]
struct CerebroSampleGroupUpdateQuery {
    // Required for access authorization in user guard middleware
    db: DatabaseId,
    project: ProjectId,
    // Sample identifier to update in this project
    id: SampleId
}

#[patch("/cerebro/samples/group")]
async fn sample_group_handler(data: web::Data<AppState>, update: web::Json<SampleGroupSchema>, query: web::Query<CerebroSampleGroupUpdateQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    if !auth_guard.user.roles.contains(&Role::Data) {
        return HttpResponse::Unauthorized().json(serde_json::json!({
            "status": "fail", "message": "You do not have permission to modify a sample", "data": serde_json::json!({"cerebro": []})
        }))
    }

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };


    match project_collection.update_many(doc! { "sample.id": &query.id }, doc! { "$set": { "sample.sample_group": &update.sample_group } }).await
    {
        Ok(_) => HttpResponse::Ok().json(serde_json::json!({
            "status": "success", "message": "Updated sample group", "data": serde_json::json!({})
        })),
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({
            "status": "error", "message": "Error updating sample group", "data": serde_json::json!({})
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let comment_schema = comment.into_inner();
    let comment = SampleComment::from_schema(&comment_schema);    
    let ids = query.id.split(",").map(|x| x.trim()).collect::<Vec<&str>>();

    match project_collection.update_many(doc! { "id":  { "$in" : &ids } }, doc! { "$push": { "sample.comments": &mongodb::bson::to_bson(&comment).unwrap() } }).await
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let ids = query.id.split(",").map(|x| x.trim()).collect::<Vec<&str>>();

    match project_collection.update_many(doc! { "id":  { "$in" : &ids } }, doc! { "$pull": { "sample.comments": { "comment_id": &query.comment_id }} }).await
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
    // Required for access authorization in middleware
    db: DatabaseId,
    project: ProjectId,
}

#[get("/cerebro/samples")]
async fn sample_handler(data: web::Data<AppState>, query: web::Query<CerebroSampleQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    match project_collection
    .distinct("sample", doc! {}) // All distinct samples in this project
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
    // Required for access authorization in middleware
    db: DatabaseId,
    project: ProjectId
}

#[get("/cerebro/workflows")]
async fn workflow_handler(data: web::Data<AppState>,query: web::Query<CerebroWorkflowQuery>, auth_guard: jwt::JwtDataMiddleware) -> HttpResponse {

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    match project_collection
    .distinct("workflow", doc! {}) // All distinct workflows in this project
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    match project_collection
    .distinct("run", doc! {}) // All distinct runs in this project
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

    let (_, _, project_collection) = match get_authorized_database_and_project_collection(&data, &query.db, &query.project, &auth_guard) {
        Ok(authorized) => authorized,
        Err(error_response) => return error_response
    };

    let pipeline = get_matched_workflow_cerebro_notaxa_pipeline(&id, &query.runs, &query.tags, &None);

    match project_collection
        .aggregate(pipeline)
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

pub fn get_authorized_database_and_project_collection(data: &web::Data<AppState>, db: &DatabaseId, project: &ProjectId, user_guard: &jwt::JwtDataMiddleware) -> Result<(MongoDatabaseName, Database, Collection<Cerebro>), HttpResponse> {
    
    let (database, project) = {
            let database_matches: Vec<&TeamDatabase> = user_guard.team.databases.iter().filter(|x| &x.id == db || &x.name == db).collect();

            if database_matches.len() != 1 {
                return Err(HttpResponse::NotFound().json(serde_json::json!({"status": "fail", "message": "Could not find requested database"})))
            }

            let database_match = database_matches[0];
            let project_matches: Vec<&ProjectCollection> = database_match.projects.iter().filter(|x| &x.id == project || &x.name == project).collect();

            if project_matches.len() != 1 {
                return Err(HttpResponse::NotFound().json(serde_json::json!({"status": "fail", "message": "Could not find requested project"})))
            }

            let project_match = project_matches[0];

            (database_match.to_owned(), project_match.to_owned())
    };

    let mongo_db = data.db.database(&database.database);
    let collection = mongo_db.collection(&project.collection);
    
    Ok((database.database, mongo_db, collection))
}

pub fn cerebro_config(cfg: &mut web::ServiceConfig, config: &Config) {
    
    cfg.service(insert_model_handler)
       .service(retrieve_model_handler) 
       .service(cerebro_identifier_handler)
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
       .service(sample_type_handler)
       .service(sample_group_handler)
       .service(history_taxa_handler)
       .service(contamination_taxa_handler_project)
    //    .service(filtered_taxa_summary_handler)
       .service(status_handler);


    if config.security.components.comments {
        // Disable comments on sample and priority taxon decision comment modification endpoints
        // Comments on priority taxon decisions and submissions are disabled in the endpoints
        cfg.service(sample_comment_handler)
           .service(delete_sample_comment_handler);
    }

    

}