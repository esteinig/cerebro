#[cfg(feature = "cli")]
#[cfg(feature = "lib")]
use typst_pdf::PdfStandards;

#[cfg(feature = "cli")]
#[cfg(feature = "lib")]
use typst::layout::Abs;

#[cfg(feature = "cli")]
#[cfg(feature = "lib")]
use std::path::Path;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;
#[cfg(target_arch = "wasm32")]
use crate::report::ReportConfig;

#[cfg(feature = "cli")]
#[cfg(feature = "lib")]
use crate::world::SystemWorld;
#[cfg(feature = "cli")]
#[cfg(feature = "lib")]
use crate::report::{PathogenDetectionReport, TrainingCompletionReport, ReportType, TemplateFormat, TemplateManager};


#[cfg(feature = "cli")]
pub struct CommandLineReportCompiler {
    world: SystemWorld,
    templates: TemplateManager,
    report_type: ReportType,
}

#[cfg(feature = "cli")]
impl CommandLineReportCompiler {
    pub fn new(root: String, report_type: ReportType) -> anyhow::Result<Self> {
        Ok(Self {
            world: SystemWorld::new(root)?,
            templates: TemplateManager::new(),
            report_type,
        })
    }
    pub fn report(&mut self, config: &Path, format: TemplateFormat) -> anyhow::Result<String>  {

        let report = match self.report_type {
            ReportType::PathogenDetection => {
                match format {
                    TemplateFormat::Json => PathogenDetectionReport::from_json(config)?,
                    TemplateFormat::Toml => PathogenDetectionReport::from_toml(config)?
                }
            },
            ReportType::TrainingCompletion => {

                match format {
                    TemplateFormat::Json => TrainingCompletionReport::from_json(config)?,
                    TemplateFormat::Toml => TrainingCompletionReport::from_toml(config)?
                }
            },
        };

        // Render the template for the selected report type
        self.templates
            .render(&self.report_type, &report).map_err(|e|{
                anyhow::anyhow!("failed to render report template: {}", e.to_string())
            })
    }
    pub fn typst(&mut self, text: &str, output: &Path) -> anyhow::Result<()> {
        
        std::fs::write(output, text).map_err(|e| anyhow::anyhow!("Failed to write Typst file: {:?}", e))
    }
    pub fn svg(&mut self, text: &str, vpath: String, output: &Path) -> anyhow::Result<()> {
        
        let document = self.world.compile(text.to_string(), vpath, None)?;

        let padding = Abs::pt(12.0);

        let svg_data = typst_svg::svg_merged(&document, padding);

        std::fs::write(output, svg_data).map_err(|e| anyhow::anyhow!("Failed to write SVG file: {:?}", e))
    }
    pub fn pdf(&mut self, text: &str, vpath: String, output: &Path) -> anyhow::Result<()> {

        let options = PdfStandards::new(&vec![
            typst_pdf::PdfStandard::V_1_7,
            typst_pdf::PdfStandard::A_2b,
        ])
        .map_err(|_| anyhow::anyhow!("failed to configure PdfStandards for Typst"))?;

        let options = typst_pdf::PdfOptions {
            ident: typst::foundations::Smart::Auto,
            timestamp: convert_datetime(chrono::Local::now()),
            page_ranges: None,
            standards: options,
        };

        {
            let document = self.world.compile(text.to_string(), vpath, None)?;

            let pdf_data = typst_pdf::pdf(&document, &options)
                .map_err(|_| anyhow::anyhow!("failed to create PDF with Typst"))?;

            std::fs::write(output, pdf_data).map_err(|e| anyhow::anyhow!("Failed to write PDF file: {:?}", e))
        }
    }
}


#[cfg(feature = "lib")]
pub struct LibraryReportCompiler {
    world: SystemWorld,
    templates: TemplateManager,
    report_type: ReportType,
}

#[cfg(feature = "lib")]
impl LibraryReportCompiler {
    pub fn new(root: String, report_type: ReportType) -> anyhow::Result<Self> {
        Ok(Self {
            world: SystemWorld::new(root)?,
            templates: TemplateManager::new(),
            report_type,
        })
    }
    pub fn report(&mut self, config: &Path, format: TemplateFormat) -> anyhow::Result<String>  {

        let report = match self.report_type {
            ReportType::PathogenDetection => {
                match format {
                    TemplateFormat::Json => PathogenDetectionReport::from_json(config)?,
                    TemplateFormat::Toml => PathogenDetectionReport::from_toml(config)?
                }
            },
            ReportType::TrainingCompletion => {

                match format {
                    TemplateFormat::Json => TrainingCompletionReport::from_json(config)?,
                    TemplateFormat::Toml => TrainingCompletionReport::from_toml(config)?
                }
            },
        };

        // Render the template for the selected report type
        self.templates
            .render(&self.report_type, &report).map_err(|e|{
                anyhow::anyhow!("failed to render report template: {}", e.to_string())
            })
    }
    pub fn typst(&mut self, text: &str, output: &Path) -> anyhow::Result<()> {
        
        std::fs::write(output, text).map_err(|e| anyhow::anyhow!("Failed to write Typst file: {:?}", e))
    }
    pub fn svg(&mut self, text: &str, vpath: String, output: &Path) -> anyhow::Result<()> {
        
        let document = self.world.compile(text.to_string(), vpath, None)?;

        let padding = Abs::pt(12.0);

        let svg_data = typst_svg::svg_merged(&document, padding);

        std::fs::write(output, svg_data).map_err(|e| anyhow::anyhow!("Failed to write SVG file: {:?}", e))
    }
    pub fn pdf(&mut self, text: &str, vpath: String, output: &Path) -> anyhow::Result<()> {

        let options = PdfStandards::new(&vec![
            typst_pdf::PdfStandard::V_1_7,
            typst_pdf::PdfStandard::A_2b,
        ])
        .map_err(|_| anyhow::anyhow!("failed to configure PdfStandards for Typst"))?;

        let options = typst_pdf::PdfOptions {
            ident: typst::foundations::Smart::Auto,
            timestamp: convert_datetime(chrono::Local::now()),
            page_ranges: None,
            standards: options,
        };

        {
            let document = self.world.compile(text.to_string(), vpath, None)?;

            let pdf_data = typst_pdf::pdf(&document, &options)
                .map_err(|_| anyhow::anyhow!("failed to create PDF with Typst"))?;

            std::fs::write(output, pdf_data).map_err(|e| anyhow::anyhow!("Failed to write PDF file: {:?}", e))
        }
    }
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub struct ReportCompiler {
    world: SystemWorld,
    templates: TemplateManager,
    report_type: ReportType,
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl ReportCompiler {
    #[wasm_bindgen(constructor)]
    pub fn new(root: String, report_type: String) -> Result<Self, JsValue> {
        Ok(Self {
            world: SystemWorld::new(root)?,
            templates: TemplateManager::new(),
            report_type: ReportType::from_str(&report_type)?,
        })
    }

    pub fn report(&mut self, config: JsValue) -> Result<String, JsValue> {

        let report = match self.report_type {
            ReportType::PathogenDetection => PathogenDetectionReport::from_js(config)?,
            ReportType::TrainingCompletion => TrainingCompletionReport::from_js(config)?
        };
        self.templates
            .render(&self.report_type, &report)
            .map_err(|e| JsValue::from(format!("failed to render report template: {}", e)))
    }

    pub fn pdf(&mut self, text: String, vpath: String) -> Result<Vec<u8>, JsValue> {
        let local_now = convert_datetime(chrono::Local::now());

        let pdf_standards = PdfStandards::new(&vec![
            typst_pdf::PdfStandard::V_1_7, 
            typst_pdf::PdfStandard::A_2b
        ]).map_err(|_| JsValue::from("failed to configure PdfStandards for Typst"))?;

        let options = typst_pdf::PdfOptions {
            ident: typst::foundations::Smart::Auto,
            timestamp: local_now,
            page_ranges: None,
            standards: pdf_standards
        };

        let pdf_data = {
            let document = self.world.compile(text, vpath, None)?;
    
            let pdf_data = typst_pdf::pdf(&document, &options)
                .map_err(|_| JsValue::from("failed to create PDF with Typst"));

            pdf_data
        };

        pdf_data
    }
    pub fn svg(&mut self, text: String, vpath: String) -> Result<Vec<String>, JsValue> {
        
        let svg_pages = {
            let document = self.world.compile(text, vpath, None)?;

            let mut pages = Vec::new();
            for page in &document.pages {
                pages.push(typst_svg::svg(page))
            }
            pages
        };

        Ok(svg_pages)
    }
}

#[cfg(feature = "cli")]
#[cfg(target_arch = "wasm32")]
/// Convert [`chrono::DateTime`] to [`typst::foundations::Datetime`]
fn convert_datetime<Tz: chrono::TimeZone>(
    date_time: chrono::DateTime<Tz>,
) -> Option<typst::foundations::Datetime> {
    typst::foundations::Datetime::from_ymd_hms(
        date_time.year(),
        date_time.month().try_into().ok()?,
        date_time.day().try_into().ok()?,
        date_time.hour().try_into().ok()?,
        date_time.minute().try_into().ok()?,
        date_time.second().try_into().ok()?,
    )
}