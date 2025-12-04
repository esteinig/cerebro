use std::{
    collections::HashSet, 
    fmt, 
    fs::File, 
    io::{BufReader, Write}, 
    path::PathBuf
};
use chrono::{NaiveDate, Utc};
use fancy_regex::Regex;
use anyhow::Result;
use serde::{Serialize,Deserialize, Deserializer};
use thiserror::Error;
use chrono::{TimeZone, SecondsFormat};

use cerebro_pipeline::{
    error::WorkflowError, 
    modules::{pathogen::PathogenDetectionTableRecord, quality::{PositiveControlConfig, PositiveControlSummary, PositiveControlSummaryBuilder, QualityControl}}, 
    nextflow::sheet::SampleSheet, 
    taxa::{filter::TaxonFilterConfig, taxon::{LineageOperations, Taxon}},
};

use crate::api::users::model::{UserId, User};
use crate::api::files::model::{SeaweedFile, SeaweedReads};

use crate::api::cerebro::schema::{
    PriorityTaxonDecisionSchema, 
    SampleCommentSchema, 
    PriorityTaxonSchema, 
    ReportSchema
};


const SCHEMA_VERSION: &str = "0.12.0";

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
    /// Represents all other cases of `WorkflowError`.
    #[error(transparent)]
    Workflow(#[from] WorkflowError),
    /// Represents all other cases of `chrono::ParseError`.
    #[error(transparent)]
    ChronoParse(#[from] chrono::ParseError),
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
    /// Indicates failure to calculate size of the model
    #[error("failed to calculate the size of the model")]
    ModelSizeCalculation,
    /// Indicates failure to construct the quality control summary from modules
    #[error("failed to construct quality control summary")]
    QualityControlSummary(#[source] WorkflowError),
}


/*
==============
Database model
==============
*/

pub type CerebroId = String;

// A struct representing a single processed workflow sample (a sequenced and analysed library)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cerebro {
    pub schema_version: String,                 // the schema version of this model
    
    pub id: CerebroId,                          // the unique identifier of this model in the database
    pub name: String,                           // the basename of the sample as processed in the workflow
    
    pub fs: Option<FileStorage>,                // the file storage configuration if running in production
    
    pub run: RunConfig,                         // the configuration of the sequence run
    pub sample: SampleConfig,                   // the configuration of the biological sample
    pub workflow: WorkflowConfig,               // the configuration of the workflow run

    pub quality: QualityControl,                // the quality control data from the parsed workflow sample

    pub taxa: Vec<Taxon>,                       // taxon data from the parsed workflow sample 
    
    // Used for internal operations, not provided to users in frontend or model downloads:

    pub tax_labels: Vec<String>                 // unique taxon labels at all taxonomic ranks detected (for searching samples in database)
    
}
impl Cerebro {
    /// Create a `Cerebro` model from a `WorkflowSample`
    /// 
    /// A model for the database represents a single sample processed through the pipeline
    /// associated with a specific biological sample, sequencing run and workflow run. 
    /// 
    /// TODO: IMPLEMENT A SIZE GUARD WHEN PARSING THE WORKFLOW SAMPLE
    /// 
    pub fn from(
        id: &str,
        quality: &QualityControl,
        taxa: &Vec<Taxon>,
        sample_sheet: Option<PathBuf>, 
        workflow_config: Option<PathBuf>,
        run_id: Option<String>,
        run_date: Option<String>,
    ) -> Result<Self, ModelError> {
            
        log::info!("Adding run and sample configs: {id}");
        let (run_config, sample_config) = match sample_sheet {
            Some(path) => {
                let sample_sheet = SampleSheet::from(&path)?;
                (RunConfig::from(&id, &sample_sheet)?, SampleConfig::from(&id, &sample_sheet)?)
            },
            None => {
                (
                    RunConfig::new(
                    &match run_id { 
                            Some(id) => id, 
                            None => String::from("Placeholder") 
                        },
                        &match run_date {
                            Some(date) => date,
                            None => Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Secs, true).to_string()
                        }
                    )?, 
                    SampleConfig::with_default(&id)?
                )
            }
        };

        log::info!("Reading workflow config: {id}");
        let workflow_config = match workflow_config {
            Some(path) => WorkflowConfig::from(&path)?,
            None => WorkflowConfig::default(),
        };

        Ok(
            Cerebro {
                schema_version: format!("{}", SCHEMA_VERSION),
                id: uuid::Uuid::new_v4().to_string(),
                name: id.to_owned(),
                fs: None,
                run: run_config, 
                sample: sample_config, 
                workflow: workflow_config,
                taxa: taxa.to_owned(),
                tax_labels: unique_taxonomic_labels(taxa),
                quality: quality.to_owned(),
            }
        )
    }
    pub fn size(&self) -> Result<usize, ModelError> {
        match serde_json::to_vec(self) {
            Ok(json_bytes) => Ok(json_bytes.len()),
            Err(_) => return Err(ModelError::ModelSizeCalculation)
        }
    }
    pub fn size_mb(&self) -> Result<f64, ModelError> {
        Ok(self.size()? as f64 / (1024.0 * 1024.0))
    }
    pub fn check_size(&self, max_size: Option<f64>, no_taxa: bool, no_taxa_max_size: Option<f64>) -> Result<(), ModelError> {

        let model_size = self.size_mb()?;
        log::info!("Model size with taxa is {:.2} MB ({})", model_size, self.name);

        if let Some(max_size) = max_size {
            if model_size > max_size {
                log::warn!("Model size ({}) >>with<< taxa is larger than allowed size ({:.2} MB)", self.name, max_size);
            }
        }

        if no_taxa {
            let mut model = self.clone();
            model.taxa.clear();
            let model_size = model.size_mb()?;
            log::info!("Model size without taxa is {:.2} MB ({})", model_size, self.name);

            if let Some(max_size) = no_taxa_max_size {
                if model_size > max_size {
                    log::warn!("Model size ({}) >>without<< taxa is larger than allowed size ({:.2} MB)", self.name, max_size);
                }
            }

        }


        Ok(())
    }
    pub fn summary_positive_control(&self) -> PositiveControlSummary {

        PositiveControlSummaryBuilder::new(
            &self.taxa, 
            &PositiveControlConfig::metagp(),
            &self.name
        ).build()
    }
    /// Pathogen detection table record per taxon
    pub fn into_pathogen_detection_table_records(&self) -> Result<Vec<PathogenDetectionTableRecord>, ModelError> {
        let mut rows = Vec::with_capacity(self.taxa.len());
        for taxon in &self.taxa {
            let row = PathogenDetectionTableRecord::from_taxon(&self.name, taxon, &taxon.evidence.alignment)?;
            rows.push(row);
        }
        Ok(rows)
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
   

}



/// Helper: Remove trailing underscore and block of uppercase letters if present.
fn get_base_name(name: &str) -> &str {
    if let Some(pos) = name.rfind('_') {
        let suffix = &name[pos + 1..];
        if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_uppercase()) {
            return &name[..pos];
        }
    }
    name
}

/// Returns a sorted vector of all unique taxonomic labels found in the input taxa.
/// For each label, also the "stripped" version (i.e. with variants removed) is added
/// this is for the 'collapse_variants' option in the filter and prevalence contamination
/// so that returned (collapsed) taxa can be searched in the hisotry data request.
pub fn unique_taxonomic_labels(taxa: &Vec<Taxon>) -> Vec<String> {
    let mut unique_labels = HashSet::new();
    
    for taxon in taxa {
        // Each label is in GTDB format like "p__Bacteria" or "s__Staphylococcus aureus"
        for label in taxon.lineage.get_labels() {
            unique_labels.insert(label.to_string());
            
            // If the label has a prefix separated by "__", split it.
            if let Some((prefix, name)) = label.split_once("__") {
                let base_name = get_base_name(name);
                let stripped_label = format!("{}__{}", prefix, base_name);
                unique_labels.insert(stripped_label);
            } else {
                // If no prefix, just add the stripped version of the whole label.
                unique_labels.insert(get_base_name(label).to_string());
            }
        }
    }
    
    let mut unique_vec: Vec<String> = unique_labels.into_iter().collect();
    unique_vec.sort();
    unique_vec
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
    pub fn from(id: &str, sample_sheet: &SampleSheet) -> Result<Self, ModelError> {

        let id = match sample_sheet.get_run_id(&id) {
            Some(run_id) => run_id,
            None => return Err(ModelError::RunIdentifier(id.to_owned()))
        };

        let date = match sample_sheet.get_run_date(&id) {
            Some(run_date) => run_date,
            None => return Err(ModelError::RunDate(id.to_owned()))
        };
        Ok(Self { id, date })
    }
    pub fn new(id: &str, date: &str) -> Result<Self, ModelError> {
        Ok(Self {
            id: id.to_string(),
            date: date.to_string()
        })
    }
}
impl Default for RunConfig {
    fn default() -> Self {
        Self {
            id: String::from("SEQUENCE-RUN"),
            date: Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Secs, true).to_string()
        }
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
    pub sample_group: Option<String>,
    pub sample_type: Option<String>,
    pub sample_date: Option<String>,
    pub comments: Vec<SampleComment>,
    pub priority: Vec<PriorityTaxon>,        // priority taxa set from the user interface
    pub reports: Vec<ReportEntry>,           // reports generated for this sample - links to a report log in the team database
    pub ercc_input_mass: Option<f64>,        // input mass in picogram
}
impl SampleConfig {
    /// Create a minimal sample configuration with defaults
    pub fn with_default(id: &str) -> Result<Self, ModelError> {

        let (mut sample_id, tags) = get_sample_regex_matches(id)?;

        log::info!("Regex parsed for id: {id}");

        if sample_id.is_empty() {
            sample_id = id.to_string()
        }

        Ok(Self {
            id: sample_id,
            tags,
            ..SampleConfig::default()
        })
    }
    /// Create a biological sample configuration from a parsed workflow sample and sample sheet
    pub fn from(id: &str, sample_sheet: &SampleSheet) -> Result<Self, ModelError> {
        let (id, tags) = get_sample_regex_matches(id)?;

        let sample_group = sample_sheet.get_sample_group(&id);

        let sample_type = sample_sheet.get_sample_type(&id);

        let sample_date = sample_sheet.get_sample_type(&id);

        let ercc_input_mass = sample_sheet.get_ercc_input(&id);

        Ok(Self{ id, tags, description: None, sample_group, sample_type, sample_date, comments: Vec::new(), priority: Vec::new(), reports: Vec::new(), ercc_input_mass })
    }
}
impl Default for SampleConfig {
    fn default() -> Self {
        Self {
            id: String::new(),
            tags: Vec::new(),
            description: None,
            sample_group: None,
            sample_type: None,
            sample_date: None,
            comments: Vec::new(),
            priority: Vec::new(),
            reports: Vec::new(),
            ercc_input_mass: None
        }
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
fn get_sample_regex_matches(file_name: &str)-> Result<(String, Vec<String>), ModelError> {
    
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
    pub workflow: String
    // pub params: WorkflowParams
}
impl WorkflowConfig {
    pub fn from(json: &PathBuf) -> Result<Self, ModelError> {
        let reader = BufReader::new(File::open(json)?);
        Ok(serde_json::from_reader(reader).map_err(ModelError::JsonSerialization)?)
    }
}
impl Default for WorkflowConfig {
    fn default() -> Self {
        Self {
            id: String::from("PLACEHOLDER"),
            name: String::from("PLACEHOLDER"),
            pipeline: String::from("PLACEHOLDER"),
            version: String::from("PLACEHOLDER"),
            started: String::from("PLACEHOLDER"),
            completed: String::from("PLACEHOLDER"),
            workflow: String::from("PLACEHOLDER")
        }
    }
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
    // pub taxon_overview: TaxonOverview,
    pub filter_config: TaxonFilterConfig,
    pub decisions: Vec<PriorityTaxonDecision>
}
impl PriorityTaxon {
    // A description of priority taxon for logging
    pub fn log_description(&self) -> String {
        format!(
            "PriorityTaxon: cerebro_ids={} id={} date={} type={} tags={} comment='{}'", 
            self.id, 
            self.date, 
            self.taxon_type, 
            // self.taxon_overview.name, 
            // self.taxon_overview.taxid, 
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
            // taxon_overview: schema.taxon_overview,
            filter_config: schema.filter_config,
            decisions: schema.decisions,
        }
    }
}


