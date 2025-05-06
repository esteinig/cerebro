use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use serde::{Deserialize, Serialize};

use cerebro_pipeline::taxa::filter::{PrevalenceContaminationConfig, TaxonFilterConfig};

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

use super::model::ModelError;

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
pub struct SampleSummarySchema {
    pub sample_ids: Vec<String>,
    pub cerebro_ids: Vec<String>,
}
impl SampleSummarySchema {
    pub fn new(sample_ids: Option<Vec<String>>) -> Self {
        Self {
            sample_ids: if let Some(sample_ids) = sample_ids { sample_ids } else { vec![] },
            cerebro_ids: vec![]
        }
    }
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

#[derive(Debug, Deserialize, Serialize)]
pub struct PatientResultSchema {
    pub review_date: String,
    pub negative: bool,
    pub organism: String,
    pub contact: String,
    pub comments: String,
    pub actions: String,
}



#[derive(Debug, Deserialize, Serialize)]
pub struct ContaminationSchema {
    pub taxid: Vec<String>,
    pub tags: Vec<String>,
    pub threshold: f64,
    pub min_rpm: f64,
    pub sample_type: Option<String>,
    pub collapse_variants: bool
}
impl ContaminationSchema {
    pub fn from_config(taxid: Vec<String>, tags: Vec<String>, config: &PrevalenceContaminationConfig) -> Self {
        Self {
            taxid,
            tags,
            threshold: config.threshold,
            min_rpm: config.min_rpm,
            sample_type: config.sample_type.clone(),
            collapse_variants: config.collapse_variants
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PrevalenceOutliers {
    pub primary: bool,
    pub secondary: bool,
    pub target: bool
}
impl PrevalenceOutliers {
    pub fn disabled() -> Self {
        Self {
            primary: false,
            secondary: false,
            target: false
        }
    }
}
impl Default for PrevalenceOutliers {
    fn default() -> Self {
        Self {
            primary: true,
            secondary: false,
            target: false
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PostFilterConfig {
    pub collapse_variants: bool,
    pub best_species: bool,
    pub best_species_min: usize,
    pub best_species_domains: Vec<String>,
    pub best_species_base_weight: Option<f64>
}
impl Default for PostFilterConfig {
    fn default() -> Self {
        Self {
            collapse_variants: true,
            best_species: true,
            best_species_min: 10,
            best_species_domains: vec![
                "Archaea".to_string(),
                "Bacteria".to_string(),
                "Eukaryota".to_string()
            ],
            best_species_base_weight: None
        }
    }
}
impl PostFilterConfig {
    pub fn disabled() -> Self {
        Self {
            collapse_variants: false,
            best_species: false,
            best_species_min: 1,
            best_species_domains: vec![],
            best_species_base_weight: None
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct MetaGpConfig {
    pub sample: String,
    pub identifiers: CerebroIdentifierSmallSchema,
    pub contamination: PrevalenceContaminationConfig,
    pub prevalence_outliers: PrevalenceOutliers,
    pub ignore_taxstr: Option<Vec<String>>
}
impl MetaGpConfig {
    pub fn new(
        sample: String, 
        controls: Option<Vec<String>>, 
        tags: Option<Vec<String>>, 
        ignore_taxstr: Option<Vec<String>>, 
        contamination: PrevalenceContaminationConfig, 
        prevalence_outliers: Option<PrevalenceOutliers>
    ) -> Self {
        Self {
            sample,
            identifiers: CerebroIdentifierSmallSchema { controls, tags },
            contamination,
            ignore_taxstr,
            prevalence_outliers: prevalence_outliers.unwrap_or(PrevalenceOutliers::default()) , // default no check for new configurations - aligns with terminal inputs
        }
    }
    pub fn with_json(
        sample: String, 
        path: &Path, 
        ignore_taxstr: Option<Vec<String>>, 
        prevalence_outliers: Option<PrevalenceOutliers>
    ) -> Result<Self, ModelError> {
        
        let rdr = BufReader::new(File::open(&path)?);
        let mut config: Self = serde_json::from_reader(rdr).map_err(|err| ModelError::JsonDeserialization(err))?;

        config.sample = sample;
        config.identifiers.assert_allowed_tags();
        
        // Only replace the values from the config file if specified

        if let Some(ignore_taxstr) = ignore_taxstr {
            config.ignore_taxstr = Some(ignore_taxstr);
        }

        if let Some(prevalence_outliers) = prevalence_outliers {
            config.prevalence_outliers = prevalence_outliers;
        }
        
        Ok(config)
    }
    pub fn to_json<P: AsRef<Path>>(&self, path: P) -> Result<(), ModelError> {
        let data = serde_json::to_string_pretty(&self).map_err(|err| ModelError::JsonSerialization(err))?;
        std::fs::write(path, data)?;
        Ok(())
    }
    pub fn from_json<P: AsRef<Path>>(path: P) -> Result<Self, ModelError> {
        let rdr = BufReader::new(File::open(&path)?);
        let config: Self = serde_json::from_reader(rdr).map_err(|err| ModelError::JsonDeserialization(err))?;
        Ok(config)
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct CerebroIdentifierSchema {
    pub sample: String,
    pub controls: Option<Vec<String>>,
    pub tags: Option<Vec<String>>
}
impl CerebroIdentifierSchema {
    pub fn new(sample: String, controls: Option<Vec<String>>, tags: Option<Vec<String>>) -> Self {
        Self {
            sample,
            controls, 
            tags
        }
    }
    pub fn from_gp_config(gp_config: &MetaGpConfig) -> Self {
        Self {
            sample: gp_config.sample.clone(),
            controls: gp_config.identifiers.controls.clone(), 
            tags: gp_config.identifiers.tags.clone()
        }
    }
    pub fn from_small_schema(sample: &str, schema: &CerebroIdentifierSmallSchema) -> Self {
        Self {
            sample: sample.to_string(),
            controls: schema.controls.clone(), 
            tags: schema.tags.clone()
        }
    }
    pub fn from_json(path: &Path) -> Result<Self, ModelError> {
        let rdr = BufReader::new(File::open(&path)?);
        serde_json::from_reader(rdr).map_err(|err| ModelError::JsonDeserialization(err))
    }
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap_or(String::new())
    }
    pub fn assert_allowed_tags(&self) {
        assert!(
            self.tags.as_ref().map_or(true, |tags| tags.iter().all(|tag| tag == "DNA" || tag == "RNA")),
            "Error: tags must only contain 'DNA' or 'RNA'"
        );
    }   
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CerebroIdentifierSmallSchema {
    pub controls: Option<Vec<String>>,
    pub tags: Option<Vec<String>>
}
impl CerebroIdentifierSmallSchema {
    pub fn new(controls: Option<Vec<String>>, tags: Option<Vec<String>>) -> Self {
        Self {
            controls, 
            tags
        }
    }
    pub fn from_json(path: &Path) -> Result<Self, ModelError> {
        let rdr = BufReader::new(File::open(&path)?);
        serde_json::from_reader(rdr).map_err(|err| ModelError::JsonDeserialization(err))
    }
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap_or(String::new())
    }
    pub fn assert_allowed_tags(&self) {
        assert!(
            self.tags.as_ref().map_or(true, |tags| tags.iter().all(|tag| tag == "DNA" || tag == "RNA")),
            "Error: tags must only contain 'DNA' or 'RNA'"
        );
    }   
}