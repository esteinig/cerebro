use serde::{Deserialize, Serialize};
use tera::{Context, Tera};

#[cfg(feature = "cli")]
use std::path::Path;
#[cfg(feature = "cli")]
use clap::ValueEnum;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::JsValue;

#[derive(Deserialize, Serialize, Clone, Debug)]
#[cfg_attr(feature = "cli", derive(ValueEnum))]
pub enum ReportType {
    PathogenDetection,
}

#[cfg(target_arch = "wasm32")]
impl ReportType {
    pub fn to_string(&self) -> String {
        match self {
            ReportType::PathogenDetection => String::from("PathogenDetection")
        }
    }
}

// Report configuration trait

pub trait ReportConfig {
    #[cfg(target_arch = "wasm32")]
    fn from_js(config: JsValue) -> Result<Self, JsValue> where Self: Sized;

    fn build_context(&self, context: &mut Context);
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
    pub contact_email: String
}


#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ReportAuthorisation {
    pub signatures: Vec<AuthorisationSignature>,
}
impl Default for ReportAuthorisation {
    fn default() -> Self {
        Self {
            signatures: vec![AuthorisationSignature::default()]
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct AuthorisationSignature {
    pub name: String,
    pub position: String,
    pub institution: String
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct ReportHeader {
    pub logo: Option<String>, // base64 encoded string
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct ReportFooter {
    pub reporting_location: String,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct ReportLegal {
    pub disclosure: String,
    pub liability: String,
    pub disclaimer: String
}

// Reports

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PathogenDetectionReport {
    pub header: ReportHeader,
    pub footer: ReportFooter,
    pub legal: ReportLegal,
    pub authorisation: ReportAuthorisation,
    pub patient_header: PatientHeader,
    pub patient_result: PatientResult
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
        }
    }
}

#[cfg(feature = "cli")]
impl PathogenDetectionReport {
    pub fn to_json_file(&self, path: &Path) -> Result<(), anyhow::Error> {
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


impl ReportConfig for PathogenDetectionReport {

    #[cfg(target_arch = "wasm32")]
    fn from_js(config: JsValue) -> Result<Self, JsValue> {
        Ok(serde_wasm_bindgen::from_value(config).map_err(|e| {
            JsValue::from(format!("Failed to deserialize config: {}", e))
        })?)
    }
    fn build_context(&self, context: &mut Context) {
        
        context.insert("report_footer_patient_id", &self.patient_header.patient_urn);
        context.insert("report_footer_date_collected", &self.patient_header.date_collected);
        context.insert("report_footer_date_reported", &self.patient_header.reporting_date);
        context.insert("report_footer_reporting_location", &self.footer.reporting_location);

        context.insert("patient_header_patient_name", &self.patient_header.patient_name);
        context.insert("patient_header_patient_dob", &self.patient_header.patient_dob);
        context.insert("patient_header_patient_urn", &self.patient_header.patient_urn);
        context.insert("patient_header_requested_doctor", &self.patient_header.requested_doctor);
        context.insert("patient_header_hospital_site", &self.patient_header.hospital_site);
        context.insert("patient_header_laboratory_number", &self.patient_header.laboratory_number);
        context.insert("patient_header_specimen_id", &self.patient_header.specimen_id);
        context.insert("patient_header_date_collected", &self.patient_header.date_collected);
        context.insert("patient_header_date_received", &self.patient_header.date_received);
        context.insert("patient_header_specimen_type", &self.patient_header.specimen_type);
        context.insert("patient_header_reporting_laboratory", &self.patient_header.reporting_laboratory);
        context.insert("patient_header_reporting_date", &self.patient_header.reporting_date);

        context.insert("patient_result_pathogen_detected", &self.patient_result.pathogen_detected);
        context.insert("patient_result_pathogen_reported", &self.patient_result.pathogen_reported);
        context.insert("patient_result_review_date", &self.patient_result.review_date);
        context.insert("patient_result_comments", &self.patient_result.comments);
        context.insert("patient_result_orthogonal_tests", &self.patient_result.orthogonal_tests);
        context.insert("patient_result_clinical_notes", &self.patient_result.clinical_notes);
        context.insert("patient_result_actions", &self.patient_result.actions);
        context.insert("patient_result_contact_name", &self.patient_result.contact_name);
        context.insert("patient_result_contact_email", &self.patient_result.contact_email);

        context.insert("report_legal_disclaimer", &self.legal.disclaimer);
        context.insert("report_legal_disclosure", &self.legal.disclosure);
        context.insert("report_legal_liability", &self.legal.liability);

        context.insert("report_authorisation_signatures", &self.authorisation.signatures);

    }
}



impl ReportType {

    #[cfg(target_arch = "wasm32")]
    // Parse the result type from a string provided through interface (WASM)
    pub fn from_str(report_type: &str) -> Result<Self, JsValue> {

        // Parse the report type
        match report_type {
            "PathogenDetection" => Ok(ReportType::PathogenDetection),
            _ => return Err(JsValue::from("invalid report type")),
        }
    }
    /// Get the name of the template for this report type
    fn template_name(&self) -> &'static str {
        match self {
            ReportType::PathogenDetection => "pathogen_detection",
        }
    }
    // Configure a report templating context
    fn context(&self, config: &dyn ReportConfig) -> Context {
        let mut context = Context::new();
        config.build_context(&mut context); // trait function configures context for each individual report type
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

        Self { tera }
    }

    /// Render a template based on the report type and context
    pub fn render(
        &self,
        report_type: &ReportType,
        report_config: &dyn ReportConfig
    ) -> Result<String, tera::Error> {
        self.tera.render(
            report_type.template_name(), 
            &report_type.context(report_config)
        )
    }
}