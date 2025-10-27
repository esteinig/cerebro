use cerebro_pipeline::{modules::{pathogen::PathogenDetectionTableRecord, quality::{QualityControl, ReadQualityControl}}, taxa::taxon::Taxon};
use serde::{Deserialize, Serialize};

use crate::api::cerebro::model::Cerebro;

use super::model::{CerebroId, RunConfig, SampleConfig, WorkflowConfig};

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CerebroIdentifierSummary {
    pub cerebro_id: String,
    pub run_id: String,
    pub run_date: String,
    pub sample_id: String,
    pub sample_tags: Vec<String>,
    pub sample_type: Option<String>,
    pub workflow_id: String,
    pub workflow_name: String,
}
#[derive(Serialize, Deserialize)]
pub struct CerebroIdentifierResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<CerebroIdentifierSummary>>
}
impl CerebroIdentifierResponse {
    pub fn success(summaries: Vec<CerebroIdentifierSummary>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Cerebro identifier summaries retrieved"),
            data: Some(summaries)
        }
    }
    pub fn sample_not_found(name: &str) -> Self {
        Self {
            status: String::from("fail"),
            message: format!("Cerebro primary sample could not be found in database: {name}"),
            data: None
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Cerebro identifier summaries could not be found in database"),
            data: None
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: None
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct FilteredTaxaResponse {
    pub status: String,
    pub message: String,
    pub data: FilteredTaxaData
}
impl FilteredTaxaResponse {
    pub fn success(taxa: Vec<Taxon>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Filtered taxa retrieved"),
            data: FilteredTaxaData { taxa }
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Filtered taxa could not be found in database"),
            data: FilteredTaxaData::default()
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: FilteredTaxaData::default()
        }
    }
}


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct RetrieveModelResponse {
    pub status: String,
    pub message: String,
    pub data: Vec<Cerebro>
}
impl RetrieveModelResponse {
    pub fn success(models: Vec<Cerebro>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Models retrieved"),
            data: models
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No models could not be found in database"),
            data: vec![]
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in model retrieval query: {}", error_message),
            data: vec![]
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PathogenDetectionTableData {
    pub records: Vec<PathogenDetectionTableRecord>,
    pub csv: String
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PathogenDetectionTableResponse {
    pub status: String,
    pub message: String,
    pub data: PathogenDetectionTableData
}
impl PathogenDetectionTableResponse {
    pub fn success(records: Vec<PathogenDetectionTableRecord>, csv: String) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Pathogen detection table data retrieved"),
            data: PathogenDetectionTableData { records, csv }
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No models could not be found in database"),
            data: PathogenDetectionTableData { records: Vec::new(), csv: String::new() }
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in pathogen detection table query: {}", error_message),
            data: PathogenDetectionTableData { records: Vec::new(), csv: String::new() }
        }
    }
}


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct QualityControlTableData {
    pub records: Vec<ReadQualityControl>,
    pub csv: String
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct QualityControlTableResponse {
    pub status: String,
    pub message: String,
    pub data: QualityControlTableData 

}
impl QualityControlTableResponse {
    pub fn success(records: Vec<ReadQualityControl>, csv: String) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Quality control table data retrieved"),
            data: QualityControlTableData { records, csv }
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("No models could not be found in database"),
            data: QualityControlTableData { records: Vec::new(), csv: String::new() }
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in quality control table query: {}", error_message),
            data: QualityControlTableData { records: Vec::new(), csv: String::new() }
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, Default)]
pub struct ContaminationTaxaData {
    pub taxid: Vec<String>
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ContaminationTaxaResponse {
    pub status: String,
    pub message: String,
    pub data: ContaminationTaxaData
}
impl ContaminationTaxaResponse {
    pub fn success(taxid: Vec<String>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Prevalence contamination taxids retrieved"),
            data: ContaminationTaxaData { taxid }
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Taxids could not be found in database"),
            data: ContaminationTaxaData::default()
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: ContaminationTaxaData::default()
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, Default)]
pub struct FilteredTaxaData {
    pub taxa: Vec<Taxon>
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct TaxonHistoryResult {
    pub id: CerebroId,
    pub sample_id: String,
    pub sample_tags: Vec<String>,
    pub sample_type: Option<String>,
    pub run_id: String,
    pub run_date: String,
    pub input_reads: u64,
    pub host_reads: Option<u64>,
    pub taxa: Vec<Taxon>,
}


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct TaxonHistoryResponse {
    pub status: String,
    pub message: String,
    pub data: Vec<TaxonHistoryResult>
}
impl TaxonHistoryResponse {
    pub fn success(data: Vec<TaxonHistoryResult>) -> Self {
        Self {
            status: String::from("success"),
            message: String::from("Taxon history retrieved"),
            data: data
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: String::from("fail"),
            message: String::from("Taxon history data could not be found in database"),
            data: vec![]
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: String::from("error"),
            message: format!("Error in database query: {}", error_message),
            data: vec![]
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SampleSummary {
    pub id: String,
    pub quality: QualityControl,
    pub run: RunConfig,
    pub sample: SampleConfig,
    pub workflow: WorkflowConfig,
    pub control_taxa: Option<Vec<Taxon>>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SampleSummaryData {
    pub summary: Vec<SampleSummary>,
    pub csv: String,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SampleSummaryResponse {
    pub status: String,
    pub message: String,
    pub data: SampleSummaryData,
}

impl SampleSummaryResponse {
    pub fn success(summary: Vec<SampleSummary>, csv: String) -> Self {
        Self {
            status: "success".to_string(),
            message: "Retrieved sample summaries for requested samples".to_string(),
            data: SampleSummaryData { summary, csv },
        }
    }

    pub fn not_found() -> Self {
        Self {
            status: "fail".to_string(),
            message: "No data found for requested samples".to_string(),
            data: SampleSummaryData {
                summary: vec![],
                csv: "".to_string(),
            },
        }
    }

    pub fn server_error(error_message: String) -> Self {
        Self {
            status: "error".to_string(),
            message: format!("Error in database query: {}", error_message),
            data: SampleSummaryData {
                summary: vec![],
                csv: "".to_string(),
            },
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CerebroResponse<T> {
    pub status: String,
    pub message: String,
    pub data: Option<T>,
}

impl<T> CerebroResponse<T> {
    pub fn ok(data: T) -> Self {
        Self { status: "ok".to_string(), message: "ok".into(), data: Some(data) }
    }
    pub fn ok_none() -> Self {
        Self { status: "ok".to_string(), message: "ok".into(), data: None }
    }
    pub fn error(message: &str) -> Self {
        Self { status: "fail".to_string(), message: message.into(), data: None }
    }
    pub fn not_found(message: &str) -> Self {
        Self { status: "fail".to_string(), message: message.into(), data: None }
    }
}