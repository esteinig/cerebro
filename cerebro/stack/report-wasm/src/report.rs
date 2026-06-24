use serde::{Deserialize, Serialize};
use tera::{Context, Tera};

#[cfg(feature = "cli")]
use clap::ValueEnum;
#[cfg(any(feature = "cli", feature = "lib"))]
use std::path::Path;
#[cfg(any(feature = "cli", feature = "lib"))]
use uuid::Uuid;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::JsValue;

#[derive(Deserialize, Serialize, Clone, Debug)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum ReportType {
    PathogenDetection,
    TrainingCompletion,
    Regression,
}

#[derive(Deserialize, Serialize, Clone, Debug)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum ReportFormat {
    Pdf,
    Svg,
    Typst,
}

#[derive(Deserialize, Serialize, Clone, Debug)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum TemplateFormat {
    Json,
    Toml,
}

#[cfg(target_arch = "wasm32")]
impl ReportType {
    pub fn to_string(&self) -> String {
        match self {
            ReportType::PathogenDetection => String::from("PathogenDetection"),
            ReportType::TrainingCompletion => String::from("TrainingCompletion"),
            ReportType::Regression => String::from("Regression"),
        }
    }
}

// Report configuration trait

pub trait ReportConfig {
    #[cfg(target_arch = "wasm32")]
    fn from_js(config: JsValue) -> Result<Self, JsValue>
    where
        Self: Sized;

    fn build_context(&self, context: &mut Context, logo_width: Option<String>);
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct PatientHeader {
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
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct PatientResult {
    pub pathogen_detected: bool,
    pub pathogen_reported: String,
    pub review_date: String,
    pub comments: String,
    pub actions: String,
    pub orthogonal_tests: String,
    pub clinical_notes: String,
    pub contact_name: String,
    pub contact_email: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ReportAuthorisation {
    pub laboratory: String,
    pub identifier: String,
    pub signatures: Vec<AuthorisationSignature>,
}
impl Default for ReportAuthorisation {
    #[cfg(target_arch = "wasm32")]
    fn default() -> Self {
        Self {
            laboratory: String::new(),
            identifier: String::new(),
            signatures: vec![AuthorisationSignature::default()],
        }
    }
    #[cfg(any(feature = "cli", feature = "lib"))]
    fn default() -> Self {
        Self {
            laboratory: String::new(),
            identifier: Uuid::new_v4().to_string(),
            signatures: vec![AuthorisationSignature::default()],
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct AuthorisationSignature {
    pub name: String,
    pub position: String,
    pub institution: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ReportHeader {
    pub logo_enabled: bool,
}
impl Default for ReportHeader {
    fn default() -> Self {
        Self { logo_enabled: true }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct ReportFooter {
    pub reporting_location: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct ReportLegal {
    pub disclosure: String,
    pub liability: String,
    pub disclaimer: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct AppendixLaboratory {
    pub enabled: bool,
    pub description: String,
    pub comments: String,
    pub header: LaboratoryHeader,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct LaboratoryHeader {
    pub protocol: String,
    pub sample_id: String,
    pub version: String,
    pub run_id: String,
    pub extraction: String,
    pub extraction_control: String,
    pub rna_depletion: String,
    pub library_control: String,
    pub adapter: String,
    pub sequencing_control: String,
    pub library: String,
    pub negative_control: String,
    pub sequencer: String,
    pub positive_control: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct AppendixBioinformatics {
    pub enabled: bool,
    pub description: String,
    pub comments: String,
    pub header: BioinformaticsHeader,
    pub libraries: String,
    pub evidence: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct BioinformaticsHeader {
    pub pipeline: String,
    pub version: String,
    pub pipeline_id: String,
    pub configuration: String,
    pub started: String,
    pub completed: String,
    pub sample_id: String,
    pub libraries: String,
    pub databases: String,
    pub taxonomy: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct BioinformaticsLibraries {}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct BioinformaticsEvidence {}

// Reports

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum MetaGptClass {
    Positive, // infectious / infectious-review
    Negative, // non-infectious / non-infectious-review
    Other,    // tumour / unknown
}
impl Default for MetaGptClass {
    fn default() -> Self {
        MetaGptClass::Other
    }
}
impl MetaGptClass {
    /// Lower-case token the template switches on for the colour-coded indicator (S8-D5).
    pub fn as_token(&self) -> &'static str {
        match self {
            MetaGptClass::Positive => "positive",
            MetaGptClass::Negative => "negative",
            MetaGptClass::Other => "other",
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct MetaGptHeader {
    pub model: String,
    pub quantization: String,
    pub parameters: String,
    pub prefetch_source: String,
    pub clinical: String,
    pub config_hash: String,
    pub sample_id: String,
    pub assay: String,
    pub call_date: String,
    pub post_filter: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct MetaGptCall {
    #[serde(default)]
    pub class: MetaGptClass,
    pub pathogen: String,
    pub statement: String,
    pub review_flagged: bool,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct MetaGptEvidenceRow {
    pub tier: String,
    pub taxon: String,
    pub taxid: String,
    pub rank: String,
    pub domain: String,
    pub rpm: String,
    pub reads: String,
    pub evidence: String,
    pub contamination: bool,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct MetaGptAssessment {
    pub enabled: bool,
    pub header: MetaGptHeader,
    pub call: MetaGptCall,
    pub evidence: Vec<MetaGptEvidenceRow>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PathogenDetectionReport {
    pub header: ReportHeader,
    pub footer: ReportFooter,
    pub legal: ReportLegal,
    pub authorisation: ReportAuthorisation,
    pub patient_header: PatientHeader,
    pub patient_result: PatientResult,
    pub appendix_laboratory: AppendixLaboratory,
    pub appendix_bioinformatics: AppendixBioinformatics,
    /// Optional generative assessment appendix (Stage 08 / S8-D1). Defaulted so existing report
    /// configs deserialize and render byte-identically when absent (S8-D6).
    #[serde(default)]
    pub meta_gpt: Option<MetaGptAssessment>,
}

impl Default for PathogenDetectionReport {
    fn default() -> Self {
        PathogenDetectionReport {
            header: ReportHeader::default(),
            footer: ReportFooter::default(),
            legal: ReportLegal::default(),
            authorisation: ReportAuthorisation::default(),
            patient_header: PatientHeader::default(),
            patient_result: PatientResult::default(),
            appendix_laboratory: AppendixLaboratory::default(),
            appendix_bioinformatics: AppendixBioinformatics::default(),
            meta_gpt: None,
        }
    }
}

#[cfg(any(feature = "cli", feature = "lib"))]
impl TrainingCompletionReport {
    pub fn to_json(&self, path: &Path) -> Result<(), anyhow::Error> {
        let json_str = serde_json::to_string_pretty(self)
            .map_err(|e| anyhow::anyhow!("failed to serialize to JSON: {}", e))?;

        std::fs::write(path, json_str)
            .map_err(|e| anyhow::anyhow!("failed to write JSON to file: {}", e))
    }
    pub fn from_json(config: &Path) -> Result<Self, anyhow::Error> {
        let file_content = std::fs::read_to_string(config)
            .map_err(|e| anyhow::anyhow!("failed to read file: {}", e))?;

        let report: Self = serde_json::from_str(&file_content)
            .map_err(|e| anyhow::anyhow!("failed to deserialize JSON: {}", e))?;

        Ok(report)
    }
    pub fn to_toml(&self, path: &Path) -> Result<(), anyhow::Error> {
        let toml_str = toml::to_string_pretty(self)
            .map_err(|e| anyhow::anyhow!("failed to serialize to TOML: {}", e))?;

        std::fs::write(path, toml_str)
            .map_err(|e| anyhow::anyhow!("failed to write TOML to file: {}", e))
    }
    pub fn from_toml(config: &Path) -> Result<Self, anyhow::Error> {
        let file_content = std::fs::read_to_string(config)
            .map_err(|e| anyhow::anyhow!("failed to read file: {}", e))?;

        let report: Self = toml::from_str(&file_content)
            .map_err(|e| anyhow::anyhow!("failed to deserialize TOML: {}", e))?;

        Ok(report)
    }
}

impl ReportConfig for PathogenDetectionReport {
    #[cfg(target_arch = "wasm32")]
    fn from_js(config: JsValue) -> Result<Self, JsValue> {
        Ok(serde_wasm_bindgen::from_value(config)
            .map_err(|e| JsValue::from(format!("Failed to deserialize config: {}", e)))?)
    }
    fn build_context(&self, context: &mut Context, logo_width: Option<String>) {
        context.insert("report_header_logo_enabled", &self.header.logo_enabled);

        if let Some(logo_width) = logo_width {
            context.insert("report_header_logo_width", &logo_width);
        } else {
            context.insert("report_header_logo_width", &String::from("11%"));
        }

        context.insert("report_footer_patient_id", &self.patient_header.patient_urn);
        context.insert(
            "report_footer_date_collected",
            &self.patient_header.date_collected,
        );
        context.insert(
            "report_footer_date_reported",
            &self.patient_header.reporting_date,
        );
        context.insert(
            "report_footer_reporting_location",
            &self.footer.reporting_location,
        );

        context.insert(
            "patient_header_patient_name",
            &self.patient_header.patient_name,
        );
        context.insert(
            "patient_header_patient_dob",
            &self.patient_header.patient_dob,
        );
        context.insert(
            "patient_header_patient_urn",
            &self.patient_header.patient_urn,
        );
        context.insert(
            "patient_header_requested_doctor",
            &self.patient_header.requested_doctor,
        );
        context.insert(
            "patient_header_hospital_site",
            &self.patient_header.hospital_site,
        );
        context.insert(
            "patient_header_laboratory_number",
            &self.patient_header.laboratory_number,
        );
        context.insert(
            "patient_header_specimen_id",
            &self.patient_header.specimen_id,
        );
        context.insert(
            "patient_header_date_collected",
            &self.patient_header.date_collected,
        );
        context.insert(
            "patient_header_date_received",
            &self.patient_header.date_received,
        );
        context.insert(
            "patient_header_specimen_type",
            &self.patient_header.specimen_type,
        );
        context.insert(
            "patient_header_reporting_laboratory",
            &self.patient_header.reporting_laboratory,
        );
        context.insert(
            "patient_header_reporting_date",
            &self.patient_header.reporting_date,
        );

        context.insert(
            "patient_result_pathogen_detected",
            &self.patient_result.pathogen_detected,
        );
        context.insert(
            "patient_result_pathogen_reported",
            &self.patient_result.pathogen_reported,
        );
        context.insert(
            "patient_result_review_date",
            &self.patient_result.review_date,
        );
        context.insert("patient_result_comments", &self.patient_result.comments);
        context.insert(
            "patient_result_orthogonal_tests",
            &self.patient_result.orthogonal_tests,
        );
        context.insert(
            "patient_result_clinical_notes",
            &self.patient_result.clinical_notes,
        );
        context.insert("patient_result_actions", &self.patient_result.actions);
        context.insert(
            "patient_result_contact_name",
            &self.patient_result.contact_name,
        );
        context.insert(
            "patient_result_contact_email",
            &self.patient_result.contact_email,
        );

        context.insert("report_legal_disclaimer", &self.legal.disclaimer);
        context.insert("report_legal_disclosure", &self.legal.disclosure);
        context.insert("report_legal_liability", &self.legal.liability);

        context.insert(
            "report_authorisation_signatures",
            &self.authorisation.signatures,
        );
        context.insert(
            "report_authorisation_laboratory",
            &self.authorisation.laboratory,
        );
        context.insert(
            "report_authorisation_identifier",
            &self.authorisation.identifier,
        );

        context.insert(
            "appendix_laboratory_enabled",
            &self.appendix_laboratory.enabled,
        );
        context.insert(
            "appendix_laboratory_description",
            &self.appendix_laboratory.description,
        );
        context.insert(
            "appendix_laboratory_comments",
            &self.appendix_laboratory.comments,
        );

        context.insert(
            "appendix_laboratory_header_protocol",
            &self.appendix_laboratory.header.protocol,
        );
        context.insert(
            "appendix_laboratory_header_sample_id",
            &self.appendix_laboratory.header.sample_id,
        );
        context.insert(
            "appendix_laboratory_header_version",
            &self.appendix_laboratory.header.version,
        );
        context.insert(
            "appendix_laboratory_header_run_id",
            &self.appendix_laboratory.header.run_id,
        );
        context.insert(
            "appendix_laboratory_header_extraction",
            &self.appendix_laboratory.header.extraction,
        );
        context.insert(
            "appendix_laboratory_header_extraction_control",
            &self.appendix_laboratory.header.extraction_control,
        );
        context.insert(
            "appendix_laboratory_header_rna_depletion",
            &self.appendix_laboratory.header.rna_depletion,
        );
        context.insert(
            "appendix_laboratory_header_library_control",
            &self.appendix_laboratory.header.library_control,
        );
        context.insert(
            "appendix_laboratory_header_adapter",
            &self.appendix_laboratory.header.adapter,
        );
        context.insert(
            "appendix_laboratory_header_sequencing_control",
            &self.appendix_laboratory.header.sequencing_control,
        );
        context.insert(
            "appendix_laboratory_header_library",
            &self.appendix_laboratory.header.library,
        );
        context.insert(
            "appendix_laboratory_header_negative_control",
            &self.appendix_laboratory.header.negative_control,
        );
        context.insert(
            "appendix_laboratory_header_sequencer",
            &self.appendix_laboratory.header.sequencer,
        );
        context.insert(
            "appendix_laboratory_header_positive_control",
            &self.appendix_laboratory.header.positive_control,
        );

        context.insert(
            "appendix_bioinformatics_enabled",
            &self.appendix_bioinformatics.enabled,
        );
        context.insert(
            "appendix_bioinformatics_description",
            &self.appendix_bioinformatics.description,
        );
        context.insert(
            "appendix_bioinformatics_comments",
            &self.appendix_bioinformatics.comments,
        );
        context.insert(
            "appendix_bioinformatics_libraries",
            &self.appendix_bioinformatics.libraries,
        );
        context.insert(
            "appendix_bioinformatics_evidence",
            &self.appendix_bioinformatics.evidence,
        );

        context.insert(
            "appendix_bioinformatics_header_pipeline",
            &self.appendix_bioinformatics.header.pipeline,
        );
        context.insert(
            "appendix_bioinformatics_header_version",
            &self.appendix_bioinformatics.header.version,
        );
        context.insert(
            "appendix_bioinformatics_header_pipeline_id",
            &self.appendix_bioinformatics.header.pipeline_id,
        );
        context.insert(
            "appendix_bioinformatics_header_configuration",
            &self.appendix_bioinformatics.header.configuration,
        );
        context.insert(
            "appendix_bioinformatics_header_started",
            &self.appendix_bioinformatics.header.started,
        );
        context.insert(
            "appendix_bioinformatics_header_completed",
            &self.appendix_bioinformatics.header.completed,
        );
        context.insert(
            "appendix_bioinformatics_header_sample_id",
            &self.appendix_bioinformatics.header.sample_id,
        );
        context.insert(
            "appendix_bioinformatics_header_libraries",
            &self.appendix_bioinformatics.header.libraries,
        );
        context.insert(
            "appendix_bioinformatics_header_databases",
            &self.appendix_bioinformatics.header.databases,
        );
        context.insert(
            "appendix_bioinformatics_header_taxonomy",
            &self.appendix_bioinformatics.header.taxonomy,
        );

        // Appendix C — generative diagnostic assessment (Stage 08 / S8-D1). Defaulted so that a
        // report without a `meta_gpt` block renders byte-identically to today (appendix disabled,
        // empty header/call/evidence).
        let gpt = self.meta_gpt.clone().unwrap_or_default();
        context.insert("appendix_metagpt_enabled", &gpt.enabled);
        context.insert("metagpt_header_model", &gpt.header.model);
        context.insert("metagpt_header_quantization", &gpt.header.quantization);
        context.insert("metagpt_header_parameters", &gpt.header.parameters);
        context.insert("metagpt_header_source", &gpt.header.prefetch_source);
        context.insert("metagpt_header_clinical", &gpt.header.clinical);
        context.insert("metagpt_header_config_hash", &gpt.header.config_hash);
        context.insert("metagpt_header_sample_id", &gpt.header.sample_id);
        context.insert("metagpt_header_assay", &gpt.header.assay);
        context.insert("metagpt_header_call_date", &gpt.header.call_date);
        context.insert("metagpt_header_post_filter", &gpt.header.post_filter);
        context.insert("metagpt_call_class", &gpt.call.class.as_token());
        context.insert("metagpt_call_pathogen", &gpt.call.pathogen);
        context.insert("metagpt_call_statement", &gpt.call.statement);
        context.insert("metagpt_call_review", &gpt.call.review_flagged);
        context.insert("metagpt_evidence", &gpt.evidence);
    }
}

// Certificate

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct TrainingCompletionReport {
    pub recipient: String,
    pub course: String,
    pub date: String,
    pub dataset: String,
    pub sensitivity: Option<f64>,
    pub specificity: Option<f64>,
}

impl Default for TrainingCompletionReport {
    fn default() -> Self {
        TrainingCompletionReport {
            recipient: String::from("Jane Doe"),
            course: String::from("Cerebro Training"),
            date: String::from("01-01-2000"),
            dataset: String::from("Cerebro Training Dataset"),
            sensitivity: None,
            specificity: None,
        }
    }
}

impl ReportConfig for TrainingCompletionReport {
    #[cfg(target_arch = "wasm32")]
    fn from_js(config: JsValue) -> Result<Self, JsValue> {
        Ok(serde_wasm_bindgen::from_value(config)
            .map_err(|e| JsValue::from(format!("Failed to deserialize config: {}", e)))?)
    }
    fn build_context(&self, context: &mut Context, logo_width: Option<String>) {
        if let Some(logo_width) = logo_width {
            context.insert("logo_width", &logo_width);
        } else {
            context.insert("logo_width", &String::from("5%"));
        }

        let sensitivity = match self.sensitivity {
            Some(v) => format!("{v:.1}%"),
            None => format!("N/A"),
        };

        let specificity = match self.specificity {
            Some(v) => format!("{v:.1}%"),
            None => format!("N/A"),
        };

        context.insert("recipient", &self.recipient);
        context.insert("course", &self.course);
        context.insert("date", &self.date);
        context.insert("dataset", &self.dataset);
        context.insert("sensitivity", &sensitivity);
        context.insert("specificity", &specificity);
    }
}

// ── Regression report (Stage 08 / S8-D3): a presentation view of the Stage-3
//    RegressionReport — populated from a stored regression report and rendered
//    by `regression_report.typ`. Mirrors how TrainingCompletion is wired. ──────

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RegressionHeader {
    pub dataset: String,
    pub dataset_version: String,
    pub baseline_id: String,
    pub model: String,
    pub quantization: String,
    pub run_id: String,
    pub evaluated: String,
    pub config_hash: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RegressionVerdict {
    pub regressed: bool,
    pub passed: bool,
    pub reasons: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RegressionStats {
    pub current_sensitivity: Option<f64>,
    pub current_specificity: Option<f64>,
    pub current_ppv: Option<f64>,
    pub current_npv: Option<f64>,
    pub baseline_sensitivity: Option<f64>,
    pub baseline_specificity: Option<f64>,
    pub baseline_ppv: Option<f64>,
    pub baseline_npv: Option<f64>,
    pub delta_sensitivity: Option<f64>,
    pub delta_specificity: Option<f64>,
    pub threshold_sensitivity: Option<f64>,
    pub threshold_specificity: Option<f64>,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RegressionMcnemar {
    pub p_value: Option<f64>,
    pub discordant: String,
    pub adjustment: String,
    pub panel_size: u64,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RegressionDecisionRow {
    pub sample_id: String,
    pub reference: String,
    pub predicted: String,
    pub decision: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RegressionReport {
    #[serde(default)]
    pub logo: bool,
    pub header: RegressionHeader,
    pub verdict: RegressionVerdict,
    pub statistics: RegressionStats,
    pub mcnemar: RegressionMcnemar,
    pub decisions: Vec<RegressionDecisionRow>,
}

/// Statistics from the evaluator are fractions in [0,1]; the report renders percentages.
fn fmt_pct(value: Option<f64>) -> String {
    match value {
        Some(v) => format!("{:.1}%", v * 100.0),
        None => String::from("N/A"),
    }
}

/// Signed percentage-point delta (current − baseline), already a fraction difference.
fn fmt_delta_pp(value: Option<f64>) -> String {
    match value {
        Some(v) => format!("{:+.1} pp", v * 100.0),
        None => String::from("N/A"),
    }
}

fn fmt_p(value: Option<f64>) -> String {
    match value {
        Some(v) => format!("{:.4}", v),
        None => String::from("N/A"),
    }
}

impl ReportConfig for RegressionReport {
    #[cfg(target_arch = "wasm32")]
    fn from_js(config: JsValue) -> Result<Self, JsValue> {
        Ok(serde_wasm_bindgen::from_value(config)
            .map_err(|e| JsValue::from(format!("Failed to deserialize config: {}", e)))?)
    }
    fn build_context(&self, context: &mut Context, logo_width: Option<String>) {
        context.insert("report_header_logo_enabled", &self.logo);
        context.insert(
            "report_header_logo_width",
            &logo_width.unwrap_or_else(|| String::from("12%")),
        );

        // Header (4-col table)
        context.insert("regression_header_dataset", &self.header.dataset);
        context.insert(
            "regression_header_dataset_version",
            &self.header.dataset_version,
        );
        context.insert("regression_header_baseline_id", &self.header.baseline_id);
        context.insert("regression_header_model", &self.header.model);
        context.insert(
            "regression_header_quantization",
            &self.header.quantization,
        );
        context.insert("regression_header_run_id", &self.header.run_id);
        context.insert("regression_header_evaluated", &self.header.evaluated);
        context.insert("regression_header_config_hash", &self.header.config_hash);

        // Verdict (colour-coded box)
        context.insert("regression_verdict_regressed", &self.verdict.regressed);
        context.insert("regression_verdict_passed", &self.verdict.passed);
        context.insert("regression_verdict_reasons", &self.verdict.reasons);

        // Statistics: current vs baseline (+ deltas, + threshold floor)
        context.insert("stats_current_sensitivity", &fmt_pct(self.statistics.current_sensitivity));
        context.insert("stats_current_specificity", &fmt_pct(self.statistics.current_specificity));
        context.insert("stats_current_ppv", &fmt_pct(self.statistics.current_ppv));
        context.insert("stats_current_npv", &fmt_pct(self.statistics.current_npv));
        context.insert("stats_baseline_sensitivity", &fmt_pct(self.statistics.baseline_sensitivity));
        context.insert("stats_baseline_specificity", &fmt_pct(self.statistics.baseline_specificity));
        context.insert("stats_baseline_ppv", &fmt_pct(self.statistics.baseline_ppv));
        context.insert("stats_baseline_npv", &fmt_pct(self.statistics.baseline_npv));
        context.insert("stats_delta_sensitivity", &fmt_delta_pp(self.statistics.delta_sensitivity));
        context.insert("stats_delta_specificity", &fmt_delta_pp(self.statistics.delta_specificity));
        context.insert("stats_threshold_sensitivity", &fmt_pct(self.statistics.threshold_sensitivity));
        context.insert("stats_threshold_specificity", &fmt_pct(self.statistics.threshold_specificity));

        // McNemar
        context.insert("mcnemar_p_value", &fmt_p(self.mcnemar.p_value));
        context.insert("mcnemar_discordant", &self.mcnemar.discordant);
        context.insert("mcnemar_adjustment", &self.mcnemar.adjustment);
        context.insert("mcnemar_panel_size", &self.mcnemar.panel_size);

        // Per-sample decisions (multi-page table via tera {% for %})
        context.insert("regression_decisions", &self.decisions);
    }
}

#[cfg(any(feature = "cli", feature = "lib"))]
impl RegressionReport {
    pub fn to_json(&self, path: &Path) -> Result<(), anyhow::Error> {
        let json_str = serde_json::to_string_pretty(self)
            .map_err(|e| anyhow::anyhow!("failed to serialize to JSON: {}", e))?;
        std::fs::write(path, json_str)
            .map_err(|e| anyhow::anyhow!("failed to write JSON to file: {}", e))
    }
    pub fn from_json(config: &Path) -> Result<Self, anyhow::Error> {
        let file_content = std::fs::read_to_string(config)
            .map_err(|e| anyhow::anyhow!("failed to read file: {}", e))?;
        let report: Self = serde_json::from_str(&file_content)
            .map_err(|e| anyhow::anyhow!("failed to deserialize JSON: {}", e))?;
        Ok(report)
    }
}

#[cfg(any(feature = "cli", feature = "lib"))]
impl PathogenDetectionReport {
    pub fn to_json(&self, path: &Path) -> Result<(), anyhow::Error> {
        let json_str = serde_json::to_string_pretty(self)
            .map_err(|e| anyhow::anyhow!("failed to serialize to JSON: {}", e))?;

        std::fs::write(path, json_str)
            .map_err(|e| anyhow::anyhow!("failed to write JSON to file: {}", e))
    }
    pub fn from_json(config: &Path) -> Result<Self, anyhow::Error> {
        let file_content = std::fs::read_to_string(config)
            .map_err(|e| anyhow::anyhow!("failed to read file: {}", e))?;

        let report: Self = serde_json::from_str(&file_content)
            .map_err(|e| anyhow::anyhow!("failed to deserialize JSON: {}", e))?;

        Ok(report)
    }
    pub fn to_toml(&self, path: &Path) -> Result<(), anyhow::Error> {
        let toml_str = toml::to_string_pretty(self)
            .map_err(|e| anyhow::anyhow!("failed to serialize to TOML: {}", e))?;

        std::fs::write(path, toml_str)
            .map_err(|e| anyhow::anyhow!("failed to write TOML to file: {}", e))
    }
    pub fn from_toml(config: &Path) -> Result<Self, anyhow::Error> {
        let file_content = std::fs::read_to_string(config)
            .map_err(|e| anyhow::anyhow!("failed to read file: {}", e))?;

        let report: Self = toml::from_str(&file_content)
            .map_err(|e| anyhow::anyhow!("failed to deserialize TOML: {}", e))?;

        Ok(report)
    }
}

impl ReportType {
    #[cfg(target_arch = "wasm32")]
    // Parse the result type from a string provided through interface (WASM)
    pub fn from_str(report_type: &str) -> Result<Self, JsValue> {
        // Parse the report type
        match report_type {
            "PathogenDetection" => Ok(ReportType::PathogenDetection),
            "TrainingCompletion" => Ok(ReportType::TrainingCompletion),
            "Regression" => Ok(ReportType::Regression),
            _ => return Err(JsValue::from("invalid report type")),
        }
    }
    /// Get the name of the template for this report type
    fn template_name(&self) -> &'static str {
        match self {
            ReportType::PathogenDetection => "pathogen_detection",
            ReportType::TrainingCompletion => "training_completion",
            ReportType::Regression => "regression_report",
        }
    }
    // Configure a report templating context
    fn context(&self, config: &dyn ReportConfig, logo_width: Option<String>) -> Context {
        let mut context = Context::new();
        config.build_context(&mut context, logo_width); // trait function configures context for each individual report type
        context
    }
}

/// Struct to manage templates
pub struct TemplateManager {
    tera: Tera,
}

impl TemplateManager {
    /// Initialize the manager with all available templates
    pub fn new() -> Self {
        let mut tera = Tera::default();

        // Add templates
        tera.add_raw_template(
            "pathogen_detection",
            include_str!("../templates/pathogen_detection.typ"),
        )
        .expect("failed to load pathogen detection template");

        tera.add_raw_template(
            "training_completion",
            include_str!("../templates/training_completion.typ"),
        )
        .expect("failed to load training completion template");

        tera.add_raw_template(
            "regression_report",
            include_str!("../templates/regression_report.typ"),
        )
        .expect("failed to load regression report template");

        Self { tera }
    }

    /// Render a template based on the report type and context
    pub fn render(
        &self,
        report_type: &ReportType,
        report_config: &dyn ReportConfig,
        logo_width: Option<String>,
    ) -> Result<String, tera::Error> {
        self.tera.render(
            report_type.template_name(),
            &report_type.context(report_config, logo_width),
        )
    }
}

#[cfg(all(test, any(feature = "cli", feature = "lib")))]
mod regression_report_tests {
    use super::*;

    fn sample_report() -> RegressionReport {
        RegressionReport {
            logo: false,
            header: RegressionHeader {
                dataset: "QC".into(),
                dataset_version: "2026.1".into(),
                baseline_id: "baseline-1".into(),
                model: "qwen3-8b-q8-0".into(),
                quantization: "q8_0".into(),
                run_id: "RUN-1".into(),
                evaluated: "2026-06-24".into(),
                config_hash: "blake3:abcd".into(),
            },
            verdict: RegressionVerdict {
                regressed: false,
                passed: true,
                reasons: "within tolerance".into(),
            },
            statistics: RegressionStats {
                current_sensitivity: Some(0.962),
                current_specificity: Some(0.948),
                current_ppv: Some(0.913),
                current_npv: Some(0.979),
                baseline_sensitivity: Some(0.958),
                baseline_specificity: Some(0.951),
                baseline_ppv: Some(0.910),
                baseline_npv: Some(0.977),
                delta_sensitivity: Some(0.004),
                delta_specificity: Some(-0.003),
                threshold_sensitivity: Some(0.90),
                threshold_specificity: Some(0.90),
            },
            mcnemar: RegressionMcnemar {
                p_value: Some(0.4807),
                discordant: "b=3, c=5".into(),
                adjustment: "continuity-corrected".into(),
                panel_size: 120,
            },
            decisions: vec![RegressionDecisionRow {
                sample_id: "BIO-001".into(),
                reference: "positive".into(),
                predicted: "infectious".into(),
                decision: "TP".into(),
            }],
        }
    }

    #[test]
    fn fmt_helpers_format_fractions_as_report_strings() {
        assert_eq!(fmt_pct(Some(0.962)), "96.2%");
        assert_eq!(fmt_pct(None), "N/A");
        assert_eq!(fmt_delta_pp(Some(0.004)), "+0.4 pp");
        assert_eq!(fmt_delta_pp(Some(-0.003)), "-0.3 pp");
        assert_eq!(fmt_p(Some(0.4807)), "0.4807");
        assert_eq!(fmt_p(None), "N/A");
    }

    #[test]
    fn build_context_emits_every_template_variable() {
        let report = sample_report();
        let mut ctx = Context::new();
        report.build_context(&mut ctx, None);

        // Every variable the regression_report.typ template references must be present.
        for key in [
            "report_header_logo_enabled",
            "report_header_logo_width",
            "regression_header_dataset",
            "regression_header_dataset_version",
            "regression_header_baseline_id",
            "regression_header_model",
            "regression_header_quantization",
            "regression_header_run_id",
            "regression_header_evaluated",
            "regression_header_config_hash",
            "regression_verdict_regressed",
            "regression_verdict_passed",
            "regression_verdict_reasons",
            "stats_current_sensitivity",
            "stats_current_specificity",
            "stats_current_ppv",
            "stats_current_npv",
            "stats_baseline_sensitivity",
            "stats_baseline_specificity",
            "stats_baseline_ppv",
            "stats_baseline_npv",
            "stats_delta_sensitivity",
            "stats_delta_specificity",
            "stats_threshold_sensitivity",
            "stats_threshold_specificity",
            "mcnemar_p_value",
            "mcnemar_discordant",
            "mcnemar_adjustment",
            "mcnemar_panel_size",
            "regression_decisions",
        ] {
            assert!(ctx.get(key).is_some(), "missing context variable: {key}");
        }

        // Spot-check a formatted value and the decisions list.
        assert_eq!(
            ctx.get("stats_current_sensitivity").unwrap().as_str().unwrap(),
            "96.2%"
        );
        assert_eq!(
            ctx.get("regression_decisions").unwrap().as_array().unwrap().len(),
            1
        );
    }

    #[test]
    fn json_round_trips() {
        let report = sample_report();
        let json = serde_json::to_string(&report).unwrap();
        let back: RegressionReport = serde_json::from_str(&json).unwrap();
        assert_eq!(back.header.model, "qwen3-8b-q8-0");
        assert_eq!(back.decisions.len(), 1);
        assert!(back.verdict.passed && !back.verdict.regressed);
    }

    #[test]
    fn report_type_round_trips_name() {
        assert_eq!(ReportType::Regression.template_name(), "regression_report");
    }

    #[test]
    fn pathogen_report_metagpt_defaults_disabled_and_back_compatible() {
        // A config without a meta_gpt block must deserialize (S8-D6) and render the appendix
        // disabled with empty header/call/evidence (byte-identical to today).
        let json = serde_json::to_string(&PathogenDetectionReport::default()).unwrap();
        let parsed: PathogenDetectionReport = serde_json::from_str(&json).unwrap();
        assert!(parsed.meta_gpt.is_none());

        let mut ctx = Context::new();
        parsed.build_context(&mut ctx, None);
        assert_eq!(
            ctx.get("appendix_metagpt_enabled").unwrap().as_bool().unwrap(),
            false
        );
        assert_eq!(ctx.get("metagpt_call_class").unwrap().as_str().unwrap(), "other");
        assert_eq!(
            ctx.get("metagpt_evidence").unwrap().as_array().unwrap().len(),
            0
        );
    }

    #[test]
    fn metagpt_class_tokens() {
        assert_eq!(MetaGptClass::Positive.as_token(), "positive");
        assert_eq!(MetaGptClass::Negative.as_token(), "negative");
        assert_eq!(MetaGptClass::Other.as_token(), "other");
    }
}
