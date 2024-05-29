use std::{
    collections::HashMap, 
    fmt, 
    io::{Write, BufReader}, 
    fs::File, 
    path::PathBuf
};
use chrono::Utc;
use fancy_regex::Regex;
use anyhow::Result;
use serde::{Serialize,Deserialize, Deserializer};
use thiserror::Error;

use cerebro_workflow::{
    error::{WorkflowError, WorkflowUtilityError}, filters::TaxonFilterConfig, module::QualityControlModule, sample::WorkflowSample, sheet::SampleSheet, taxon::{Taxon, TaxonOverview}
};

use crate::api::users::model::{UserId, User};
use crate::api::files::model::{SeaweedFile, SeaweedReads};

use crate::api::cerebro::schema::{
    PriorityTaxonDecisionSchema, 
    SampleCommentSchema, 
    PriorityTaxonSchema, 
    ReportSchema
};


const SCHEMA_VERSION: &str= "0.8.0";

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum ModelError {
    /// Indicates failure with JSON serialization.
    #[error("failed to serialize JSON")]
    JsonSerialization(#[source] serde_json::Error),
    /// Indicates failure with JSON serialization.
    #[error("failed to deserialize JSON")]
    JsonDeserialization(#[source] serde_json::Error),
    /// Represents all other cases of `std::io::Error`.
    #[error(transparent)]
    IOError(#[from] std::io::Error),
    /// Represents all other cases of `workflow::utils::WorkflowUtilityError`.
    #[error(transparent)]
    WorkflowUtility(#[from] WorkflowUtilityError),
    /// Represents failure to find a run identifier in the sample sheet
    #[error("failed to find run identifier in sample sheet for sample identifier: {0}")]
    RunIdentifier(String),
    /// Represents failure to find a run date in the sample sheet
    #[error("failed to find run date in sample sheet for sample identifier: {0}")]
    RunDate(String),
    /// Represents failure to find a run date in the sample sheet
    #[error("failed to find sample group in sample sheet for sample identifier: {0}")]
    SampleGroup(String),
    /// Represents failure to find a sample type in the sample sheet
    #[error("failed to find sample type in sample sheet for sample identifier: {0}")]
    SampleType(String),
    /// Indicates failure to compile the sample id regex
    #[error("failed to compile sample identifier regex")]
    SampleIdRegex(#[source] fancy_regex::Error),
    /// Indicates failure to obtain the sample id regex capture
    #[error("failed to get sample identifier regex capture")]
    SampleIdRegexCapture(#[source] fancy_regex::Error),
    /// Indicates failure to obtain the sample id regex capture match
    #[error("failed to get sample identifier regex capture match")]
    SampleIdRegexCaptureMatch,
    /// Indicates failure to compile the sample tag regex
    #[error("failed to compile sample tag regex")]
    SampleTagRegex(#[source] fancy_regex::Error),
    /// Indicates failure to obtain the sample tag regex capture
    #[error("failed to get sample tag regex capture")]
    SampleTagRegexCapture(#[source] fancy_regex::Error),
    /// Indicates failure to obtain the sample tag regex capture match
    #[error("failed to get sample tag regex capture match")]
    SampleTagRegexCaptureMatch,
    /// Indicates failure to construct the quality control summary from modules
    #[error("failed to construct quality control summary")]
    QualityControlSummary(#[source] WorkflowError),
}


/*
========================
Run configuration
========================
*/

pub type RunId = String;

// A struct representing a sequence run configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunConfig {
    pub id: RunId,
    pub date: String,
}
impl RunConfig {
    pub fn from(sample_sheet: &SampleSheet, workflow_sample: &WorkflowSample) -> Result<Self, ModelError> {

        let id = match sample_sheet.get_run_id(&workflow_sample.id) {
            Some(run_id) => run_id,
            None => return Err(ModelError::RunIdentifier(workflow_sample.id.to_owned()))
        };

        let date = match sample_sheet.get_run_date(&workflow_sample.id) {
            Some(run_date) => run_date,
            None => return Err(ModelError::RunDate(workflow_sample.id.to_owned()))
        };
        Ok(Self { id, date })
    }
}

/*
========================
Sample configuration
========================
*/

pub type SampleId = String;
pub type Tag = String;

// A struct representing a biological sample data configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleConfig {
    pub id: SampleId,                       // identifier of the sample i.e. file base name stripped of tags
    pub tags: Vec<Tag>,
    pub description: Option<String>,
    pub sample_group: String,
    pub sample_type: String,
    pub comments: Vec<SampleComment>,
    pub priority: Vec<PriorityTaxon>,        // priority taxa set from the user interface
    pub reports: Vec<ReportEntry>,           // reports generated for this sample - links to a report log in the team database
    pub ercc_input_mass: Option<f64>,        // input mass in picogram
}
impl SampleConfig {
    /// Create a biological sample configuration from a parsed workflow sample and sample sheet
    pub fn from(sample_sheet: &SampleSheet, workflow_sample: &WorkflowSample) -> Result<Self, ModelError> {
        let (id, tags) = get_sample_regex_matches(&workflow_sample.id)?;

        let sample_group = match sample_sheet.get_sample_group(&workflow_sample.id) {
            Some(sample_group) => sample_group,
            None => return Err(ModelError::SampleGroup(workflow_sample.id.to_owned()))
        };

        let sample_type = match sample_sheet.get_sample_type(&workflow_sample.id) {
            Some(sample_type) => sample_type,
            None => return Err(ModelError::SampleGroup(workflow_sample.id.to_owned()))
        };

        let ercc_input_mass = sample_sheet.get_ercc_input(&workflow_sample.id);

        Ok(Self{ id, tags, description: None, sample_group, sample_type, comments: Vec::new(), priority: Vec::new(), reports: Vec::new(), ercc_input_mass })
    }
     /// Create a biological sample configuration from a parsed workflow sample
     pub fn from_sample(workflow_sample: &WorkflowSample) -> Result<Self, ModelError> {
        let (id, tags) = get_sample_regex_matches(&workflow_sample.id)?;
        Ok(Self{ id, tags, description: None, sample_group: String::new(), sample_type: String::new(), comments: Vec::new(), priority: Vec::new(), reports: Vec::new(), ercc_input_mass: None })
    }
}

pub type ReportId = String;

// A struct representing an entry for a report to be stored
// in the models that were reported on
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportEntry {
    pub id: ReportId,
    pub date: String,
    pub user_id: UserId,
    pub user_name: String,
    pub negative: bool,
    pub organism: String,
    pub review_date: String,
    pub report_text: Option<String>,
    pub report_pdf: Option<bool>
}
impl ReportEntry {
    pub fn from_schema(id: String, schema: &ReportSchema, report_text: Option<String>, is_pdf: Option<bool>) -> Self {
        Self {
            id,
            date: Utc::now().to_string(),
            user_id: schema.user_id.clone(),
            user_name: schema.user_name.clone(),
            negative: schema.patient_result.negative,
            organism: schema.patient_result.organism.clone(),
            review_date: schema.patient_result.review_date.clone(),
            report_text,
            report_pdf: is_pdf
        }
    }
    pub fn log_description(&self) -> String {
        format!(
            "Report entry added: report_id={} report_date={} user_id={} user_name='{}' organism='{}' review_date='{}' negative={}", 
            self.id, self.date, self.user_id, self.user_name, self.organism, self.review_date, self.negative
        )
    }
}

// A struct representing a biological sample data configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleComment {
    pub user_id: UserId,
    pub user_name: String,
    pub user_positions: Vec<String>,
    pub user_title: Option<String>,
    pub comment_id: String,
    pub comment_date: String,
    pub comment_text: String
}
impl SampleComment {
    pub fn from_schema(comment: &SampleCommentSchema) -> Self {
        Self {
            user_id: comment.user_id.clone(),
            user_name: comment.user_name.clone(),
            user_positions: comment.user_positions.clone(),
            user_title: comment.user_title.clone(),
            comment_id: uuid::Uuid::new_v4().to_string(),
            comment_date: Utc::now().to_string(),
            comment_text: comment.comment_text.clone()
        }
    }
    pub fn log_description(&self) -> String {
        format!("Sample comment added: comment_id={} comment_date={} comment_text='{}'", self.comment_id, self.comment_date, self.comment_text)
    }
}

// Utility function to extract the biological sample identifier and library tags [strict]
fn get_sample_regex_matches(file_name: &String)-> Result<(String, Vec<String>), ModelError> {
    
    let mut sample_id = String::new();
    let sample_id_regex = Regex::new("^([^_]+)(?=__)").map_err(ModelError::SampleIdRegex)?;
    for caps in sample_id_regex.captures_iter(&file_name) {
        let sample_name = caps.map_err(ModelError::SampleIdRegexCapture)?.get(0);

        sample_id = match sample_name {
            Some(name) => name.as_str().to_string(),
            None => return Err(ModelError::SampleIdRegexCaptureMatch)
        };
        break; // always use the first capture
    }

    let mut library_tags: Vec<String> = Vec::new();
    let sample_id_regex = Regex::new("(?<=__)([A-Za-z0-9]+)").map_err(ModelError::SampleTagRegex)?;
    for caps in sample_id_regex.captures_iter(&file_name) {
        let tags = caps.map_err(ModelError::SampleTagRegexCapture)?.get(0);
        let tag = match tags {
            Some(name) => name.as_str().to_string(),
            None => return Err(ModelError::SampleTagRegexCaptureMatch)
        };
        library_tags.push(tag)
    }
    Ok((sample_id, library_tags))
}




// A struct representing the file system access to 
// the input and workflow files for this library
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileStorage {
    pub raw_reads: Option<SeaweedReads>,          // input reads including  host background
    pub qc_reads: Option<SeaweedReads>,           // reads after qc and background depletion - add additional human depletion after contig assembly
    pub wf_data: Option<SeaweedFile>,             // for now compressed archive - implement more differentiated file get + contig/classified read access late
}


/*
========================
Workflow configuration
========================
*/

pub type WorkflowId = String;

// Output format from workflow defined in: lib/utils.nf
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowConfig {
    pub id: WorkflowId,
    pub name: String,
    pub pipeline: String,
    pub version: String,
    pub started: String,
    pub completed: String,
    pub params: WorkflowParams
}
impl WorkflowConfig {
    pub fn from(json: &PathBuf) -> Result<Self, ModelError> {
        let reader = BufReader::new(File::open(json)?);
        Ok(serde_json::from_reader(reader).map_err(ModelError::JsonSerialization)?)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParams {
    pub production: WorkflowParamsProduction,
    pub qc: WorkflowParamsQc
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQc {
    pub enabled: bool,
    pub deduplication: WorkflowParamsQcDeduplication,
    pub reads: WorkflowParamsQcReads,
    pub controls: WorkflowParamsQcControls,
    pub host: WorkflowParamsQcHost
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsProduction {
    pub enabled: bool,
    pub sample_sheet: PathBuf,
    pub api: WorkflowParamsProductionApi,
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsProductionApi {
    pub enabled: bool,
    pub url: String,
    pub token: String,
    pub upload: WorkflowParamsProductionApiUpload
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsProductionApiUpload {
    pub enabled: bool,
    pub team: Option<String>,
    pub db: Option<String>,
    pub project: Option<String>,
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum DeduplicationMethod {
    #[serde(alias="none")]
    None, 
    #[serde(alias="naive")]
    Naive,
    #[serde(alias="umi-naive")]
    UmiNaive,
    #[serde(alias="umi-calib")]
    UmiCalib,
    #[serde(alias="umi-tools")]
    UmiTools,
    #[serde(alias="umi-tools-naive")]
    UmiToolsNaive
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcDeduplication {
    pub enabled: bool,
    pub method: DeduplicationMethod,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcReads {
    pub fastp: WorkflowParamsQcFastp,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcFastp {
    pub enabled: bool,
    pub min_read_length: Option<u64>,
    pub cut_tail_quality: Option<u64>,
    pub complexity_threshold: Option<u64>,
    pub adapter_auto_detect: bool,
    pub adapter_file: Option<PathBuf>,
    pub adapter_seq_1: Option<String>,
    pub adapter_seq_2: Option<String>,
    pub trim_poly_g: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcControls {
    pub ercc: WorkflowParamsQcErcc,
    pub phage: WorkflowParamsQcPhage,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcErcc {
    pub enabled: bool,
    pub fasta: Option<PathBuf>
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcPhage {
    pub enabled: bool,
    pub fasta: Option<PathBuf>,
    pub identifiers: WorkflowParamsQcPhageIdentifiers,
}



#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcPhageIdentifiers {
    pub dna_extraction: Option<String>,
    pub rna_extraction: Option<String>,
    pub sequencing: Option<String>,
}



#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcHost {
    pub depletion: WorkflowParamsQcHostDepletion
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct WorkflowParamsQcHostDepletion {
    pub enabled: bool,
    pub databases: String,
    pub references: String,
    pub taxa: String,
    pub direct: String,
    pub min_cov: u64,
    pub min_len: u64,
    pub min_mapq: u64
}

// Actix Web fails when specifying custom serde deserializers - it might be because
// they are used to deserialize from JSON into struct, then the struct is serialized
// into JSON for the request, and the request is again deserialized into a struct at
// the endpoint and failing because the value is not serialized back into its original
// format when adding the struct JSON as request payload.

fn _split_comma_path_str<'de, D>(deserializer: D) -> Result<Vec<PathBuf>, D::Error>
where
    D: Deserializer<'de>,
{
    let str_sequence = String::deserialize(deserializer)?;
    Ok(str_sequence
        .split(',')
        .map(|x| PathBuf::from(x.trim()))
        .collect())
}

fn _split_whitespace_str<'de, D>(deserializer: D) -> Result<Vec<String>, D::Error>
where
    D: Deserializer<'de>,
{
    let str_sequence = String::deserialize(deserializer)?;
    Ok(str_sequence
        .split_whitespace()
        .map(|item| item.trim().to_owned())
        .collect())
}



/*
========================
User configuration
========================
*/

pub type PriorityTaxonId = String;
pub type PriorityTaxonDecisionId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TaxonType {
    Unknown,
    Contaminant,
    Pathogen
}
impl fmt::Display for TaxonType {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            TaxonType::Unknown => write!(f, "Unknown"),
            TaxonType::Contaminant => write!(f, "Contaminant"),
            TaxonType::Pathogen => write!(f, "Pathogen")
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DecisionType {
    Accept,
    Reject
}
impl fmt::Display for DecisionType {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DecisionType::Accept => write!(f, "accept"),
            DecisionType::Reject => write!(f, "reject"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PriorityTaxonDecision {
    pub id: PriorityTaxonDecisionId,
    pub date: String,
    pub user_id: UserId,
    pub user_name: String,
    pub decision: DecisionType,
    pub comment: String,
}
impl PriorityTaxonDecision {
    pub fn from(schema: &PriorityTaxonDecisionSchema, user: &User) -> Self {
        Self {
            id: uuid::Uuid::new_v4().to_string(),
            date: Utc::now().to_string(),
            user_id: user.id.to_owned(),
            user_name: user.name.to_owned(),
            decision: schema.decision.to_owned(),
            comment: schema.decision_comment.to_owned()
        }
    }
    pub fn log_description(&self, schema: &PriorityTaxonDecisionSchema) -> String {
        format!(
            "PriorityTaxonDecision: priority_taxon_id={} id={} date={} decision={} type={} name='{}' taxid={} comment='{}'", 
            schema.id, self.id, self.date, self.decision, schema.taxon_type, schema.taxon_name, schema.taxon_taxid, schema.decision_comment
        )
    }
}

// User configurations set in the front-end interface
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PriorityTaxon {
    pub id: PriorityTaxonId,
    pub user_name: String,
    pub user_id: UserId,
    pub comment: String,
    pub date: String,
    pub evidence_tags: Vec<String>,            // joined tags, not arrays
    pub cerebro_identifiers: Vec<CerebroId>,   
    pub taxon_type: TaxonType,
    pub taxon_overview: TaxonOverview,
    pub filter_config: TaxonFilterConfig,
    pub decisions: Vec<PriorityTaxonDecision>
}
impl PriorityTaxon {
    // A description of priority taxon for logging
    pub fn log_description(&self) -> String {
        format!(
            "PriorityTaxon: cerebro_ids={} id={} date={} type={} name='{}' taxid={} tags={} comment='{}'", 
            self.id, 
            self.date, 
            self.taxon_type, 
            self.taxon_overview.name, 
            self.taxon_overview.taxid, 
            self.evidence_tags.join(","),
            self.cerebro_identifiers.join(","),
            self.comment
        )
    }
    // Create from schema
    pub fn from_schema(schema: PriorityTaxonSchema) -> Self {
        Self {
            id: uuid::Uuid::new_v4().to_string(),
            date: Utc::now().to_string(),
            user_id: schema.user_id,
            user_name: schema.user_name,
            comment: schema.comment,
            evidence_tags: schema.evidence_tags,          // joined tags, not arrays
            cerebro_identifiers:schema.cerebro_identifiers,  
            taxon_type: schema.taxon_type,
            taxon_overview: schema.taxon_overview,
            filter_config: schema.filter_config,
            decisions: schema.decisions,
        }
    }
}


/*
==============
Database model
==============
*/

// quality_summary: Vec::from([QualityControlRow::from(&workflow_sample.qc_module)])

pub type CerebroId = String;

// A struct representing a single processed workflow sample (a sequenced and analysed library)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cerebro {
    pub schema_version: String,             // the schema version of this model

    pub id: CerebroId,                      // the unique identifier of this model in the database
    pub name: String,                       // the basename of the sample as processed in the workflow

    pub fs: Option<FileStorage>,            // the file storage configuration if running in production

    pub run: RunConfig,                     // the configuration of the sequence run
    pub sample: SampleConfig,               // the configuration of the biological sample
    pub workflow: WorkflowConfig,           // the configuration of the workflow run

    pub taxa: HashMap<String, Taxon>,               // the dictionary of taxonomic identifiers and taxon data from the parsed workflow sample (legacy dictionary, could be simple list)
    pub quality: QualityControlModule,              // the quality control data from the parsed workflow sample
}
impl Cerebro {
    /// Create a `Cerebro` model from a `WorkflowSample`
    /// 
    /// A model for the database represents a single sample processed through the pipeline
    /// associated with a specific biological sample, sequencing run and workflow run. 
    /// 
    /// IMPLEMENT A SIZE GUARD WHEN PARSING THE WORKFLOW SAMPLE
    /// 
    pub fn from_workflow_sample(
        workflow_sample: &WorkflowSample, 
        sample_sheet: &PathBuf, 
        workflow_config: &PathBuf
    ) -> Result<Self, ModelError> {
            Ok(
                Cerebro {
                    schema_version: format!("{}", SCHEMA_VERSION),
                    id: uuid::Uuid::new_v4().to_string(),
                    name: workflow_sample.id.to_owned(),
                    fs: None,
                    run:  RunConfig::from(&SampleSheet::from(&sample_sheet)?, &workflow_sample)?, 
                    sample: SampleConfig::from(&SampleSheet::from(&sample_sheet)?, &workflow_sample)?, 
                    workflow: WorkflowConfig::from(workflow_config)?,
                    taxa: workflow_sample.taxa.to_owned(),
                    quality: workflow_sample.qc_module.to_owned(),
                }
            )
    }
    pub fn from_json(file: &PathBuf) -> Result<Self, ModelError>{
        let model: Cerebro = serde_json::from_reader(File::open(&file)?).map_err(ModelError::JsonDeserialization)?;
        Ok(model)
    }
    pub fn write_json(&self, file: &PathBuf) -> Result<(), ModelError>{
        let mut file = File::create(&file)?;
        let json_string = serde_json::to_string_pretty(&self).map_err(ModelError::JsonSerialization)?;
        write!(file, "{}", json_string)?;
        Ok(())
    }
    pub fn update_sample_id(&self, sample_id: &str, sample_tags: Option<Vec<String>>) -> Self {
       
       // TODO: Updates sample identifier throughout model
        
        log::warn!("Updating sample identifier for model: {} ({})", self.name, self.id);

        let mut cc = self.clone();

        cc.name = sample_id.to_string();
        cc.sample.id = sample_id.to_string();

        if let Some(tags) = sample_tags {
            cc.sample.tags = tags;
        }
        
        cc.quality.id = sample_id.to_string();

        cc.taxa = self.taxa.iter().map(|(taxid, taxon)| {
            let taxon = taxon.update_evidence_sample_id(sample_id);
            (taxid.to_owned(), taxon)
        }).collect::<HashMap<String, Taxon>>();

        cc
    }

}



