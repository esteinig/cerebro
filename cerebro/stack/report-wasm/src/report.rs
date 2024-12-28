use serde::Deserialize;
use serde_wasm_bindgen::from_value;
use tera::{Context, Tera};
use wasm_bindgen::JsValue;

/// Enum representing available report templates
pub enum ReportType {
    PathogenDetection,
}

// Report configuration for trait

pub trait ReportConfig {
    fn from_js(config: JsValue) -> Result<Self, JsValue> where Self: Sized;
    fn build_context(&self, context: &mut Context);
}

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

pub struct PatientResult {
    pub pathogen_detected: bool,
    pub pathogen_reported: String,
    pub review_date: String,
    pub comments: String,
    pub actions: String,
    pub contact_name: String,
    pub contact_email: String
}


#[derive(Deserialize)]
pub struct ReportAuthorisation {
    pub review_date: String,
    pub signatures: Vec<AuthorisationSignature>,
}

#[derive(Clone, Deserialize)]
pub struct AuthorisationSignature {
    pub name: String,
    pub position: String,
    pub institution: String

}

#[derive(Deserialize)]
pub struct ReportHeader {
    pub logo: Option<String>, // base64 encoded string
}

#[derive(Deserialize)]
pub struct ReportFooter {
    pub patient_id: String,
    pub date_collected: String,
    pub date_reported: String,
    pub reporting_location: String,
}

#[derive(Deserialize)]
pub struct ReportLegal {
    pub disclosure: String,
    pub liability: String,
    pub disclaimer: String
}
impl ReportLegal {
    pub fn new(disclosure: &str, disclaimer: &str, liability: &str) -> Self {
        Self {
            disclosure: disclosure.to_string(),
            disclaimer: disclaimer.to_string(),
            liability: liability.to_string()
        }
    }
}

// Reports for ReportType


#[derive(Deserialize)]
pub struct PathogenDetectionReport {
    pub header: ReportHeader,
    pub footer: ReportFooter,
    pub legal: ReportLegal,
    pub authorisation: ReportAuthorisation
}
impl ReportConfig for PathogenDetectionReport {
    fn from_js(config: JsValue) -> Result<Self, JsValue> {
        Ok(from_value(config).map_err(|e| {
            JsValue::from(format!("Failed to deserialize config: {}", e))
        })?)
    }
    fn build_context(&self, context: &mut Context) {
        
    }
}



impl ReportType {
    // Parse the result type from a string provided through interface (WASM)
    pub fn from_str(report_type: &str) -> Result<Self, JsValue> {

        // Parse the report type
        match report_type {
            "PathogenDetection" => Ok(ReportType::PathogenDetection),
            _ => return Err(JsValue::from("Invalid report type")),
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
            include_str!("../templates/pathogen_detection.tera"),
        )
        .expect("Failed to load pathogen detection template");

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