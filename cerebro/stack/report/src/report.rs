
use serde::{Serialize, Deserialize};
use handlebars::Handlebars;
use std::path::PathBuf;
use std::fs;

use typst::foundations::Smart;
use typst::{eval::Tracer, layout::Abs};

use cerebro_pipeline::quality::QualityControlSummary;
use cerebro_model::{
    api::cerebro::model::{CerebroId, PriorityTaxonDecision}, 
    api::cerebro::schema::ReportSchema, 
};

use crate::typst::TypstWrapperWorld;
use crate::error::ReportError;

// Template sub-sections - can be read from TOML

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct  ReportPatient {
    toml: String,
    header: PatientHeader,
    result: PatientResult,

    #[serde(skip_deserializing)]
    pub toml_path: Option<PathBuf>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PatientHeader {
    pub toml: String,

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

    #[serde(skip_deserializing)]
    pub toml_path: Option<PathBuf>
}
impl PatientHeader {
    pub fn from_schema(schema: &ReportSchema) -> Self {
        Self {
            toml: String::new(),
            toml_path: None,
            patient_name: schema.patient_header.patient_name.clone(),
            patient_urn: schema.patient_header.patient_urn.clone(),
            patient_dob: schema.patient_header.patient_dob.clone(),
            requested_doctor: schema.patient_header.requested_doctor.clone(),
            hospital_site: schema.patient_header.hospital_site.clone(),
            laboratory_number: schema.patient_header.laboratory_number.clone(),
            specimen_id: schema.patient_header.specimen_id.clone(),
            date_collected: schema.patient_header.date_collected.clone(),
            date_received: schema.patient_header.date_received.clone(),
            specimen_type: schema.patient_header.specimen_type.clone(),
            reporting_laboratory: schema.patient_header.reporting_laboratory.clone(),
            reporting_date: schema.patient_header.reporting_date.clone(),
            reporting_location: schema.patient_header.reporting_location.clone(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PatientResult {
    pub toml: String,

    pub review_date: String,
    pub negative: bool,
    pub organism: String,
    pub contact: String,
    pub comments: String,
    pub actions: String,

    #[serde(skip_deserializing)]
    pub toml_path: Option<PathBuf>
}
impl PatientResult {
    pub fn from_schema(schema: &ReportSchema) -> Self {
        Self {
            toml: String::new(),
            toml_path: None,
            review_date: schema.patient_result.review_date.clone(),
            negative: schema.patient_result.negative.clone(),
            organism: schema.patient_result.organism.clone(),
            contact: schema.patient_result.contact.clone(),
            comments: schema.patient_result.comments.clone(),
            actions: schema.patient_result.actions.clone(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TemplateConfig {
    pub file: String,
    pub logo: String,
    pub name: String,
    pub version: String,
    pub description: String,
    pub legal: TemplateLegal,
    pub authorisation: TemplateAuthorisation,

    #[serde(skip_deserializing)]
    pub file_path: PathBuf,
    #[serde(skip_deserializing)]
    pub logo_path: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TemplateLegal {
    pub disclosure: String,
    pub liability: String,
    pub disclaimer: String
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TemplateAuthorisation {
    pub issuing_laboratory_text: String,
    pub report_identifier_text: String,
    pub report_identifier: Option<String>,
    pub signature: Vec<AuhtorisationSignature>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AuhtorisationSignature {
    pub name: String,
    pub position: String,
    pub institution: String
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportSample {
    pub toml: String,

    pub sample_id: String,
    pub run_id: String,
    pub negative: bool,
    pub laboratory: SampleLaboratory,
    pub bioinformatics: SampleBioinformatics,
    pub audit: SampleAudit,

    // Legal disclaimer so it's available in the rendering loop, 
    // weirdly the template engine does not handle the outer 
    // scope variables inside the loop
    #[serde(skip_deserializing)]
    pub disclaimer: String,

    #[serde(skip_deserializing)]
    pub toml_path: Option<PathBuf>
}
impl ReportSample {
    pub fn from_api(
        schema: &ReportSchema,
        template: &TemplateConfig, 
        assay_template: &AssayTemplate, 
        bioinformatics_template: &BioinformaticsTemplate,
        quality_summaries: &Vec<QualityControlSummary>
    ) -> Self {
        
        Self {
            toml: String::new(),
            toml_path: None,
            disclaimer: template.legal.disclaimer.clone(),

            sample_id: schema.sample_id.to_string(),
            run_id: schema.run_id.to_string(),
            negative: schema.patient_result.negative,

            laboratory: SampleLaboratory { 
                toml: String::new(),
                toml_path: None,
                comments: schema.laboratory_comments.clone(),
                description: assay_template.description.clone(),
                header: assay_template.header.clone()
            },
            bioinformatics: SampleBioinformatics { 
                toml: String::new(),
                toml_path: None,

                description: bioinformatics_template.description.clone(), 
                comments: schema.bioinformatics_comments.clone(), 
                header: BioinformaticsHeader { 
                    pipeline: schema.workflow.pipeline.clone(), 
                    version: schema.workflow.version.clone(), 
                    run_id: schema.workflow.id.clone(), 
                    started: schema.workflow.started.clone(), 
                    completed: schema.workflow.completed.clone(), 
                    libraries: quality_summaries.iter().map(
                        |summary| match &summary.library_tag { Some(tag) => tag.to_string(), None => String::from("ERROR") }
                    ).collect::<Vec<String>>().join(", "), 
                    configuration: String::from("CNS-ASSAY"), 
                    databases: String::from("CIPHER 0.7.0"), 
                    taxonomy: String::from("NCBI")
                }, 
                library: quality_summaries.iter().map(
                    |summary| BioinformaticsLibrary::from_summary(summary)
                ).collect(), 
                evidence: BioinformaticsEvidence::from_schema(&schema)
            },
            audit: SampleAudit { 
                include: false, 
                decisions: Vec::new(), 
                models: Vec::new()
            }
        }
    }
}



#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssayTemplate {
    pub description: String,
    pub header: LaboratoryHeader
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BioinformaticsTemplate {
    pub description: String
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleLaboratory {
    pub toml: String,
    pub description: String,
    pub comments: String,
    pub header: LaboratoryHeader,

    #[serde(skip_deserializing)]
    pub toml_path: Option<PathBuf>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LaboratoryHeader {
    pub protocol: String,
    pub version: String,
    pub extraction: String,
    pub ribosomal: String,
    pub adapters: String,
    pub library: String,
    pub sequencer: String,
    pub extraction_control: String,
    pub library_control: String,
    pub sequencing_control: String,
    pub negative_control: String,
    pub positive_control: String
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleBioinformatics {
    pub toml: String,
    pub description: String,
    pub comments: String,
    pub header: BioinformaticsHeader,
    pub library: Vec<BioinformaticsLibrary>,
    pub evidence: Vec<BioinformaticsEvidence>,

    #[serde(skip_deserializing)]
    pub toml_path: Option<PathBuf>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BioinformaticsHeader {
    pub pipeline: String,
    pub version: String,
    pub run_id: String,
    pub started: String,
    pub completed: String,
    pub libraries: String,
    pub configuration: String,
    pub databases: String,
    pub taxonomy: String
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BioinformaticsLibrary {
    pub tag: String,
    pub reads: String,
    pub reads_qc: String,
    pub extraction_control: String,
    pub library_control: String,
    pub sequencing_control: String
}
impl BioinformaticsLibrary {
    pub fn from_summary(qc_summary: &QualityControlSummary) -> Self {
        Self {
            tag: match &qc_summary.library_tag { Some(tag) => tag.to_string(), None => "ERROR".to_string() },
            reads: qc_summary.total_reads.to_string(),
            reads_qc: qc_summary.output_reads.to_string(),
            extraction_control: match &qc_summary.control_status_dna_extraction { Some(status) => status.to_string(), None => String::from("PASS") },
            sequencing_control: match &qc_summary.control_status_sequencing { Some(status) => status.to_string(), None => String::from("N/A") },
            library_control: match &qc_summary.control_status_library { Some(status) => status.to_string(), None => String::from("PASS") }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BioinformaticsEvidence {
    pub organism: String,
    pub taxid: String,
    pub rpm: String,
    pub contigs: String,
    pub negative_control: bool,
    pub other: Option<String>,
}
impl BioinformaticsEvidence {
    pub fn from_schema(schema: &ReportSchema) -> Vec<Self> {
        
        match &schema.priority_taxon {
            Some(priority_taxon) => {
                vec![Self {
                    organism: priority_taxon.taxon_overview.name.clone(),
                    taxid: priority_taxon.taxon_overview.taxid.clone(),
                    rpm: format!("{:.1}", priority_taxon.taxon_overview.rpm),
                    contigs: priority_taxon.taxon_overview.contigs.to_string(),
                    negative_control: match schema.priority_taxon_negative_control { Some(v) => v, None => false },
                    other: schema.priority_taxon_other_evidence.clone()
                }]
            },
            None => Vec::new()
        }
        
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleAudit {
    include: bool,
    decisions: Vec<PriorityTaxonDecision>,
    models: Vec<CerebroId>
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportConfig {
    pub template: TemplateConfig,
    pub patient: ReportPatient,
    pub sample: Vec<ReportSample>
}   


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClinicalReport {
    pub id: uuid::Uuid,
    pub config: ReportConfig
}
impl ClinicalReport {
    pub fn from_toml(base_config: &PathBuf, patient_config: Option<PathBuf>, patient_header_config: Option<PathBuf>, patient_result_config: Option<PathBuf>, sample_configs: Option<Vec<PathBuf>>) -> Result<Self, ReportError>  {
        
        let toml_str = std::fs::read_to_string(base_config).map_err(|err| ReportError::TomlConfigFile(err))?;
        let config: ReportConfig = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;

        let report_id = uuid::Uuid::new_v4();
        let mut report_config = Self { id: report_id, config };

        // Adjusts relative paths from input configuration
        report_config.config.template.file_path = get_adjusted_path(&report_config.config.template.file, base_config);
        report_config.config.template.logo_path = get_adjusted_path(&report_config.config.template.logo, base_config);

        // Patient constructors

        match (patient_config, report_config.config.patient.toml.as_str()) {
            (None, "") => {
                log::info!("No separate patient config file has been provided");
                report_config.config.patient.toml_path = None;
            },
            (None, _) =>  {
                log::info!("A separate patient config file has been provided in the base template");
                let path = get_adjusted_path(&report_config.config.patient.toml, base_config);
                let toml_str = std::fs::read_to_string(path.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                let config: ReportPatient = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                
                report_config.config.patient = config;
                report_config.config.patient.toml_path = Some(path);
            },
            (Some(patient_file), _) => {

                log::info!("A separate patient config file has been provided by the method execution");
                let toml_str = std::fs::read_to_string(patient_file.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                let config: ReportPatient = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                
                report_config.config.patient = config;
                report_config.config.patient.toml_path = Some(patient_file);
            }
        };

        match (patient_header_config, report_config.config.patient.header.toml.as_str()) {
            (None, "") => {
                log::info!("No separate patient header config file has been provided");
                report_config.config.patient.header.toml_path = None;
            },
            (None, _) =>  {
                log::info!("A separate patient header config file has been provided in the base template");
                let path = get_adjusted_path(&report_config.config.patient.header.toml, base_config);
                let toml_str = std::fs::read_to_string(path.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                let config: PatientHeader = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                
                report_config.config.patient.header = config;
                report_config.config.patient.header.toml_path = Some(path);
            },
            (Some(patient_header_file), _) => {

                log::info!("A separate patient header config file has been provided by the method execution");
                let toml_str = std::fs::read_to_string(patient_header_file.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                let config: PatientHeader = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                
                report_config.config.patient.header = config;
                report_config.config.patient.header.toml_path = Some(patient_header_file);
            }
        };

        match (patient_result_config, report_config.config.patient.result.toml.as_str()) {
            (None, "") => {
                log::info!("No separate patient header config file has been provided");
                report_config.config.patient.result.toml_path = None;
            },
            (None, _) =>  {
                log::info!("A separate patient header config file has been provided in the base template");
                let path = get_adjusted_path(&report_config.config.patient.result.toml, base_config);
                let toml_str = std::fs::read_to_string(path.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                let config: PatientResult = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                
                report_config.config.patient.result = config;
                report_config.config.patient.result.toml_path = Some(path);
            },
            (Some(patient_result_file), _) => {

                log::info!("A separate patient header config file has been provided by the method execution");
                let toml_str = std::fs::read_to_string(patient_result_file.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                let config: PatientResult = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                
                report_config.config.patient.result = config;
                report_config.config.patient.result.toml_path = Some(patient_result_file);
            }
        };

        // Sample constructors

        match sample_configs {
            Some(sample_config_files) => {
                log::info!("Separate sample config files have been provided by the method execution");
                let mut samples = Vec::new();
                for path in sample_config_files {
                    let toml_str = std::fs::read_to_string(path.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                    let mut config: ReportSample = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                    
                    config.toml_path = Some(path);
                    samples.push(config)

                }
                report_config.config.sample = samples;

            },
            None => {
                let mut samples = Vec::new();
                for mut sample in report_config.config.sample {
                    match sample.toml.as_str() {
                        "" => {
                            log::info!("No separate sample config file has been provided for sample {}", sample.sample_id);
                            sample.toml_path = None;
                        },
                        _ => {
                            log::info!("A separate sample config file has been provided in the base template execution");
                            let path = get_adjusted_path(&sample.toml, base_config);
                            let toml_str = std::fs::read_to_string(path.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                            let config: ReportSample = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                            sample = config;
                            sample.toml_path = Some(path);
                        }
                    }
                    samples.push(sample)
                };
                report_config.config.sample = samples;
            }
        }

        // Nested laboratory and bioinformatics section if provided as files in toml field paths

        let mut samples = Vec::new();
        for mut sample_config in report_config.config.sample.clone() {
            
            match sample_config.laboratory.toml.as_str() {
                "" => {
                    log::info!("No separate laboratory config file has been provided for {}", sample_config.sample_id);
                    sample_config.laboratory.toml_path = None;
                },
                _ =>  {
                    log::info!("A separate laboratory config file has been provided for {}",sample_config.sample_id);
                    let path = get_adjusted_path(&sample_config.laboratory.toml, base_config);
                    let toml_str = std::fs::read_to_string(path.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                    let config: SampleLaboratory = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                    
                    sample_config.laboratory = config;
                    sample_config.laboratory.toml_path = Some(path);
                },
                
            };

            match sample_config.bioinformatics.toml.as_str() {
                "" => {
                    log::info!("No separate bioinformatics config file has been provided for {}", sample_config.sample_id);
                    sample_config.bioinformatics.toml_path = None;
                },
                _ =>  {
                    log::info!("A separate bioinformatics config file has been provided for {}", sample_config.sample_id);
                    let path = get_adjusted_path(&sample_config.bioinformatics.toml, base_config);
                    let toml_str = std::fs::read_to_string(path.clone()).map_err(|err| ReportError::TomlConfigFile(err))?;
                    let config: SampleBioinformatics = toml::from_str(&toml_str).map_err(|err| ReportError::TomlConfigParse(err))?;
                    
                    sample_config.bioinformatics = config;
                    sample_config.bioinformatics.toml_path = Some(path);
                },  
            };

            sample_config.disclaimer = report_config.config.template.legal.disclaimer.clone();

            samples.push(sample_config);
        }

        report_config.config.sample = samples;

        // Set the ClinicalReport UUID as the templated report identifier if the option is not present in the parsed template
        if let None = report_config.config.template.authorisation.report_identifier {
            report_config.config.template.authorisation.report_identifier = Some(report_id.to_string());
        }
        


        Ok(report_config)
    }
    // Configures the clinical report in API
    pub fn from_api(
        schema: &ReportSchema, 
        template: &mut TemplateConfig, 
        assay_template: &AssayTemplate, 
        bioinformatics_template: &BioinformaticsTemplate, 
        quality_summaries: &Vec<QualityControlSummary> 
    ) -> Self {
        

        // Always set the ClinicalReport UUID as the report identifier via API
        let report_id = uuid::Uuid::new_v4();
        template.authorisation.report_identifier = Some(report_id.to_string());

        Self {
            id: report_id,
            config: ReportConfig { 
                template: template.clone(),
                patient: ReportPatient { 
                    toml: String::new(), 
                    header: PatientHeader::from_schema(&schema), 
                    result: PatientResult::from_schema(&schema), 
                    toml_path: None
                },
                sample: vec![
                    ReportSample::from_api(schema, template, assay_template, bioinformatics_template, quality_summaries)
                ]
            }
        }
    }
    pub fn render_pdf(&self, output: Option<PathBuf>, output_template: Option<PathBuf>, output_svg: Option<PathBuf>) -> Result<Vec<u8>, handlebars::RenderError> {

        let current_dir: PathBuf = std::env::current_dir()?;

        let template_string = self.render_template(output_template)?;

        // Create Typst world with template as content:
        let world = TypstWrapperWorld::new(current_dir.display().to_string(), template_string);

        // Render document
        let mut tracer = Tracer::default();
        let document = typst::compile(&world, &mut tracer).expect("Error compiling typst.");

        let pdf = typst_pdf::pdf(&document, Smart::Auto, None);

        if let Some(path) = output {
                fs::write(&path, &pdf).expect("Error writing PDF.");
                
        }
        if let Some(path) = output_svg {
            let svg_string = typst_svg::svg_merged(&document, Abs::pt(2.0));
            fs::write(&path, &svg_string).expect("Error writing SVG.");
        }

        Ok(pdf)
    }

    pub fn render_template(&self, output: Option<PathBuf>) -> Result<String, handlebars::RenderError> {

        let mut handlebars = Handlebars::new();

        handlebars.register_template_file(&self.config.template.name, &self.config.template.file_path)?;
        let content_template = handlebars.render(&self.config.template.name, &self.config)?;

        if let Some(file) = output {
            fs::write(file.clone(), content_template.clone()).expect(
                &format!("Unable to write template to file: {}", file.display())
            );
        }

        Ok(content_template)

    } 
}


pub fn get_adjusted_path(path: &String, config_file: &PathBuf) -> PathBuf {
    match path.starts_with("/") {
        // Interpret as absolute path
        true => PathBuf::from(path.to_owned()),
        // Interpret as relative path to the main report configuration file (TOML)
        false => {
            let mut full_path = config_file.parent().expect(
                &format!("Failed to get parent directory of report configuration file {}", &config_file.display())
            ).to_path_buf();
            full_path.push(path);
            full_path
        }
    }
}