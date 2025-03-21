use serde::{Deserialize, Serialize};

use cerebro_pipeline::taxa::filter::TaxonFilterConfig;

use crate::api::users::model::UserId;
use crate::api::cerebro::model::{
    PriorityTaxonId, 
    PriorityTaxon, 
    DecisionType, 
    TaxonType, 
    CerebroId, 
    PriorityTaxonDecision, 
    WorkflowConfig
};

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PriorityTaxonDecisionSchema {
    pub id: PriorityTaxonId,
    pub decision: DecisionType,
    pub decision_comment: String,
    pub taxon_name: String,
    pub taxon_taxid: String,
    pub taxon_type: TaxonType
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ReportData {
}


#[derive(Debug, Deserialize)]
pub struct SampleCommentSchema {
    pub user_id: UserId,
    pub user_name: String,
    pub user_positions: Vec<String>,
    pub user_title: Option<String>,
    pub comment_text: String,
}


#[derive(Deserialize)]
pub struct SampleDeleteSchema {
    pub sample_id: Vec<String>,
}


#[derive(Deserialize, Serialize, Debug)]
pub struct SampleSummaryQcSchema {
    pub sample_ids: Vec<String>,
    pub cerebro_ids: Vec<String>,
}



#[derive(Deserialize)]
pub struct SampleDescriptionSchema {
    pub description: String
}

#[derive(Deserialize)]
pub struct SampleTypeSchema {
    pub sample_type: String
}

#[derive(Deserialize)]
pub struct SampleGroupSchema {
    pub sample_group: String
}



#[derive(Deserialize, Serialize)]
pub struct TaxaSummarySchema {
    pub sample_ids: Vec<String>,
    pub run_ids: Vec<String>,
    pub workflow_ids: Vec<String>,
    pub workflow_names: Vec<String>,
    pub filter_config: TaxonFilterConfig
}

#[derive(Deserialize, Serialize)]
pub struct PriorityTaxonSchema {
    pub user_name: String,
    pub user_id: UserId,
    pub comment: String,
    pub evidence_tags: Vec<String>,            // joined tags, not arrays
    pub cerebro_identifiers: Vec<CerebroId>,   
    pub taxon_type: TaxonType,
    // pub taxon_overview: TaxonOverview,
    pub filter_config: TaxonFilterConfig,
    pub decisions: Vec<PriorityTaxonDecision>
}

#[derive(Deserialize, Serialize)]
pub struct ReportSchema {
    pub ids: Vec<CerebroId>,
    pub user_id: UserId,
    pub user_name: String,
    pub sample_id: String,
    pub run_id: String,
    pub workflow: WorkflowConfig,
    pub patient_header: PatientHeaderSchema,
    pub patient_result: PatientResultSchema,
    pub laboratory_comments: String,
    pub bioinformatics_comments: String,
    pub priority_taxon: Option<PriorityTaxon>,
    pub priority_taxon_negative_control: Option<bool>,
    pub priority_taxon_other_evidence: Option<String>,
}

#[derive(Deserialize, Serialize)]
pub struct PatientHeaderSchema {
    pub patient_name: String,
    pub patient_urn: String,
    pub patient_dob: String,
    pub requested_doctor: String,
    pub hospital_site: String,
    pub laboratory_number: String,
    pub specimen_id: String,
    pub date_collected: String,
    pub date_received: String,
    pub specimen_type: String,
    pub reporting_laboratory: String,
    pub reporting_date: String,
    pub reporting_location: String,
}

#[derive(Deserialize, Serialize)]
pub struct PatientResultSchema {
    pub review_date: String,
    pub negative: bool,
    pub organism: String,
    pub contact: String,
    pub comments: String,
    pub actions: String,
}



#[derive(Deserialize)]
pub struct ContaminationSchema {
    pub taxid: Vec<String>,
    pub tags: Vec<String>,
    pub threshold: f64,
    pub min_rpm: f64,
    pub sample_type: Option<String>,
}

