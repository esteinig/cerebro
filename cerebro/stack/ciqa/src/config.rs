use std::{fs::File, io::BufReader, path::PathBuf};

use cerebro_pipeline::taxa::filter::{PrevalenceContaminationConfig, TaxonFilterConfig};
use serde::{Deserialize, Serialize};

use crate::{error::CiqaError, terminal::EvaluateArgs};


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EvaluationConfig {
    pub controls: Vec<String>,
    pub tags: Vec<String>,
    pub filter: TaxonFilterConfig,
    pub prevalence: PrevalenceContaminationConfig
}
impl EvaluationConfig {
    pub fn from_json(path: &PathBuf) -> Result<Self, CiqaError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let schema: Self = serde_json::from_reader(reader)?;
        schema.assert_allowed_tags();
        Ok(schema)
    }
    pub fn from_evaluate_args(args: &EvaluateArgs) -> Self {
        let schema = Self {
            controls: args.controls.clone(),
            tags: args.tags.clone(),
            filter: TaxonFilterConfig::validation(),
            prevalence: PrevalenceContaminationConfig::validation(),
        };
        schema.assert_allowed_tags();
        schema
    }
    pub fn assert_allowed_tags(&self) {
        assert!(
            self.tags.iter().all(|tag| tag == "DNA" || tag == "RNA"),
            "Error: tags must only contain 'DNA' or 'RNA'"
        );
    }   
}
