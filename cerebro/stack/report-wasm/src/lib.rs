use chrono::{Datelike, Timelike};
use typst_pdf::PdfStandards;
use wasm_bindgen::prelude::*;

mod world;

use crate::world::SystemWorld;

#[wasm_bindgen]
pub struct ReportCompiler {
    world: SystemWorld,
}

#[wasm_bindgen]
impl ReportCompiler {
    #[wasm_bindgen(constructor)]
    pub fn new(root: String, request_data: &js_sys::Function) -> Self {

        Self {
            world: SystemWorld::new(root, request_data)
        }
    }

    pub fn pdf(&mut self, text: String, path: String) -> Result<Vec<u8>, JsValue> {
                
        let local_now = convert_datetime(chrono::Local::now());
        let pdf_standards = PdfStandards::new(&vec![
            typst_pdf::PdfStandard::V_1_7, 
            typst_pdf::PdfStandard::A_2b
        ]).map_err(|_| JsValue::from("Failed to configure PdfStandards for Typst"))?;

        let options = typst_pdf::PdfOptions {
            ident: typst::foundations::Smart::Auto,
            timestamp: local_now,
            page_ranges: None,
            standards: pdf_standards
        };

        let document = self
            .world
            .compile(text, path)
            .map_err(|_| JsValue::from("Failed to compile document with Typst"))?;

        typst_pdf::pdf(&document, &options)
            .map_err(|_| JsValue::from("Failed to create PDF from compiled document with Typst"))
            
    }

    pub fn add_font(&mut self, data: Vec<u8>) {
        self.world.add_font(data);
    }
}

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