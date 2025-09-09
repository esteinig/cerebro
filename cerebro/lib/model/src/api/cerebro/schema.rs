use std::collections::HashSet;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use serde::{Deserialize, Serialize};
use std::io::{Write, BufWriter};
use cerebro_pipeline::taxa::{taxon::Taxon};
use cerebro_pipeline::taxa::filter::{TargetList, TaxonFilterConfig};

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
pub struct QualityControlTableSchema {
    pub sample_ids: Vec<String>,
    pub cerebro_ids: Vec<String>,
}
impl QualityControlTableSchema {
    pub fn new(sample_ids: Vec<String>) -> Self {
        Self {
            sample_ids,
            cerebro_ids: vec![]
        }
    }
}

#[derive(Deserialize, Serialize, Debug)]
pub struct PathogenDetectionTableSchema {
    pub sample_ids: Vec<String>,
    pub cerebro_ids: Vec<String>
}
impl PathogenDetectionTableSchema {
    pub fn new(sample_ids: Vec<String>) -> Self {
        Self {
            sample_ids,
            cerebro_ids: vec![]
        }
    }
}
impl Default for PathogenDetectionTableSchema {
    fn default() -> Self {
        PathogenDetectionTableSchema { sample_ids: vec![], cerebro_ids: vec![] }
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



#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrevalenceContaminationConfig {
    pub threshold: f64,
    pub min_rpm: f64,
    pub sample_type: Option<String>,
    pub collapse_variants: bool,
    pub outliers: PrevalenceOutliers
}
impl Default for PrevalenceContaminationConfig {
    fn default() -> Self {
        Self {
            threshold: 0.60,
            min_rpm: 0.0,
            sample_type: None,
            collapse_variants: false,
            outliers: PrevalenceOutliers::default()
        }
    }
}
impl PrevalenceContaminationConfig {
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
    pub fn none() -> Self {
        Self {
            threshold: 0.0,
            min_rpm: 0.0,
            sample_type: None,
            collapse_variants: false,
            outliers: PrevalenceOutliers::disabled()
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
    pub best_species_base_weight: Option<f64>,
    pub exclude_phage: bool,
    pub exclude_phage_list: HashSet<String>
}
impl Default for PostFilterConfig {
    fn default() -> Self {
        Self {
            collapse_variants: true,
            best_species: true,
            best_species_min: 3,
            best_species_domains: vec![
                "Archaea".to_string(),
                "Bacteria".to_string(),
                "Eukaryota".to_string()
            ],
            best_species_base_weight: None,
            exclude_phage: true,
            exclude_phage_list: HashSet::from_iter(
                TargetList::gp_prokaryote_viruses().to_vec()
            )
        }
    }
}
impl PostFilterConfig { 
    pub fn with_default(collapse_variants: bool, best_species: bool, min_species: usize, species_domains: Vec<String>, exclude_phage: bool) -> Self {
        let mut config = Self::default();
        config.collapse_variants = collapse_variants;
        config.best_species = best_species;
        config.best_species_min = min_species;
        config.best_species_domains = species_domains;
        config.exclude_phage = exclude_phage;
        config
    }
    pub fn disabled() -> Self {
        Self {
            collapse_variants: false,
            best_species: false,
            best_species_min: 1,
            best_species_domains: vec![],
            best_species_base_weight: None,
            exclude_phage: false,
            exclude_phage_list: HashSet::new()
        }
    }
}


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct TieredFilterConfig {
    pub primary: TaxonFilterConfig,
    pub secondary: TaxonFilterConfig,
    pub target: TaxonFilterConfig
}
impl TieredFilterConfig {
    pub fn with_ignore_taxstr(&self, ignore_taxstr: Vec<String>) -> Self {
        let mut config = self.clone();
        config.primary.ignore_taxstr = match config.primary.ignore_taxstr {
            Some(v) => Some([v, ignore_taxstr.clone()].concat().to_vec()),
            None => Some(ignore_taxstr.clone())
        };
        config.secondary.ignore_taxstr = match config.secondary.ignore_taxstr {
            Some(v) => Some([v, ignore_taxstr.clone()].concat().to_vec()),
            None => Some(ignore_taxstr.clone())
        };
        config.target.ignore_taxstr = match config.target.ignore_taxstr {
            Some(v) => Some([v, ignore_taxstr.clone()].concat().to_vec()),
            None => Some(ignore_taxstr.clone())
        };
        config
    }
    pub fn default(ignore_taxstr: Option<Vec<String>>) -> Self {
        Self {
            primary: TaxonFilterConfig::gp_above_threshold(
                ignore_taxstr.clone()
            ),
            secondary: TaxonFilterConfig::gp_below_threshold(
                ignore_taxstr.clone()
            ),
            target: TaxonFilterConfig::gp_target_threshold(
                ignore_taxstr.clone()
            )
        }
    }
    pub fn none() -> Self {
        Self {
            primary: TaxonFilterConfig::default(),
            secondary: TaxonFilterConfig::default(),
            target: TaxonFilterConfig::default()
        }
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

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct MetaGpConfig {
    pub sample: String,
    pub identifiers: CerebroIdentifierSmallSchema,
    pub filter_configs: TieredFilterConfig,
    pub contamination: PrevalenceContaminationConfig,
}
impl MetaGpConfig {
    pub fn new(
        sample: String, 
        controls: Option<Vec<String>>, 
        tags: Option<Vec<String>>, 
        filter_configs: TieredFilterConfig,
        contamination: PrevalenceContaminationConfig, 
    ) -> Self {
        Self {
            sample,
            identifiers: CerebroIdentifierSmallSchema { controls, tags },
            contamination,
            filter_configs,
        }
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


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PrefetchData {
    pub primary: Vec<Taxon>,
    pub secondary: Vec<Taxon>,
    pub target: Vec<Taxon>,
    pub primary_contamination: Vec<Taxon>,
    pub secondary_contamination: Vec<Taxon>,
    pub target_contamination: Vec<Taxon>,
    pub primary_filter: TaxonFilterConfig,
    pub secondary_filter: TaxonFilterConfig,
    pub target_filter: TaxonFilterConfig,
    pub config: MetaGpConfig
}
impl PrefetchData {
    pub fn new(
        primary: Vec<Taxon>,
        secondary: Vec<Taxon>,
        target: Vec<Taxon>,
        primary_contamination: Vec<Taxon>,
        secondary_contamination: Vec<Taxon>,
        target_contamination: Vec<Taxon>,
        primary_filter: &TaxonFilterConfig,
        secondary_filter: &TaxonFilterConfig,
        target_filter: &TaxonFilterConfig,
        config: &MetaGpConfig
    ) -> Self {
        Self {
            primary,
            secondary,
            target,
            primary_contamination,
            secondary_contamination,
            target_contamination,
            primary_filter: primary_filter.clone(),
            secondary_filter: secondary_filter.clone(),
            target_filter: target_filter.clone(),
            config: config.clone()
        }
    }
    /// Remove from `secondary` any Taxon whose `lineage` also appears in `primary`,
    /// then remove from `target` any Taxon whose `lineage` also appears in the 
    /// pruned `secondary` and 'primary' - including prevalence contamination
    pub fn prune(&mut self) {
        
        // collect all lineages present in primary
        let primary_lineages: HashSet<_> =
            self.primary.iter().map(|t| t.lineage.clone()).collect();
        let primary_contam_lineages: HashSet<_> =
            self.primary_contamination.iter().map(|t| t.lineage.clone()).collect();

        // drop from secondary anything already in primary
        self.secondary
            .retain(|t| !primary_lineages.contains(&t.lineage));
        self.secondary
            .retain(|t| !primary_contam_lineages.contains(&t.lineage));

        // now collect lineages in the (pruned) secondary
        let secondary_lineages: HashSet<_> =
            self.secondary.iter().map(|t| t.lineage.clone()).collect();
        let secondary_contam_lineages: HashSet<_> =
            self.secondary_contamination.iter().map(|t| t.lineage.clone()).collect();

        // drop from target anything already in secondary
        self.target
            .retain(|t| !secondary_lineages.contains(&t.lineage));
        self.target
            .retain(|t| !secondary_contam_lineages.contains(&t.lineage));

        // drop from target anything already in primary
        self.target
            .retain(|t| !primary_lineages.contains(&t.lineage));
        self.target
            .retain(|t| !primary_contam_lineages.contains(&t.lineage));

    }
    pub fn to_json(&self, path: &Path) -> Result<(), ModelError> {
        let data = serde_json::to_string_pretty(self).map_err(ModelError::JsonSerialization)?;
        let mut writer = BufWriter::new(File::create(path)?);
        write!(writer, "{data}")?;
        Ok(())
    }
    pub fn from_json<P: AsRef<Path>>(path: P) -> Result<Self, ModelError> {
        let data = std::fs::read_to_string(path)?;
        let result = serde_json::from_str::<PrefetchData>(&data).map_err(ModelError::JsonDeserialization)?;
        Ok(result)
    }
}