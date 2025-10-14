use std::{fs::File, io::{BufWriter, Write}, path::Path};
use cerebro_model::api::{cerebro::schema::TestResult, training::model::TrainingResultRecord};
use env_logger::Builder;
use env_logger::fmt::Color;
use log::{LevelFilter, Level};
use serde::{Deserialize, Serialize};

use crate::error::HttpClientError;


pub fn init_logger() {

    Builder::new()
        .format(|buf, record| {
            let timestamp = buf.timestamp();

            let mut red_style = buf.style();
            red_style.set_color(Color::Red).set_bold(true);
            let mut green_style = buf.style();
            green_style.set_color(Color::Green).set_bold(true);
            let mut white_style = buf.style();
            white_style.set_color(Color::White).set_bold(false);
            let mut orange_style = buf.style();
            orange_style.set_color(Color::Rgb(255, 102, 0)).set_bold(true);
            let mut apricot_style = buf.style();
            apricot_style.set_color(Color::Rgb(255, 195, 0)).set_bold(true);

            let msg = match record.level(){
                Level::Warn => (orange_style.value(record.level()), orange_style.value(record.args())),
                Level::Info => (green_style.value(record.level()), white_style.value(record.args())),
                Level::Debug => (apricot_style.value(record.level()), apricot_style.value(record.args())),
                Level::Error => (red_style.value(record.level()), red_style.value(record.args())),
                _ => (white_style.value(record.level()), white_style.value(record.args()))
            };

            writeln!(
                buf,
                "{} [{}] - {}",
                white_style.value(timestamp),
                msg.0,
                msg.1
            )
        })
        .filter(None, LevelFilter::Info)
        .init();
}

// Matches the META-GPT implementation - need to find a better way to
// provide this struct across Cerebro 


#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub enum Diagnosis {
    Infectious,
    InfectiousReview,
    NonInfectious,
    NonInfectiousReview,
    Tumor,
    Unknown
}

#[derive(Debug, Serialize, Deserialize)]
pub struct DiagnosticResult {
    pub diagnosis: Diagnosis,
    pub candidates: Vec<String>,
    pub pathogen: Option<String>,
}
impl DiagnosticResult {
    pub fn non_infectious() -> Self {
        Self {
            diagnosis: Diagnosis::NonInfectious,
            candidates: vec![],
            pathogen: None
        }
    }
    pub fn to_json(&self, path: &Path) -> Result<(), HttpClientError> {
        let agent_state = serde_json::to_string_pretty(self).map_err(|err| HttpClientError::SerdeFailure(err))?;
        let mut writer = BufWriter::new(File::create(path)?);
        write!(writer, "{agent_state}")?;
        Ok(())
    }
    pub fn from_json<P: AsRef<Path>>(path: P) -> Result<Self, HttpClientError> {
        let data = std::fs::read_to_string(path)?;
        let result = serde_json::from_str::<DiagnosticResult>(&data)?;
        Ok(result)
    }

    pub fn from_training_result(record: &TrainingResultRecord) -> Self {
        
        let candidates = match &record.candidates {
            Some(candidates) => candidates.split(";").map(|x| x.trim().to_string()).collect(),
            None => vec![]
        };
        
        Self {
            diagnosis: match record.result { 
                TestResult::Positive => Diagnosis::Infectious,
                TestResult::Negative => Diagnosis::NonInfectious
            },
            pathogen: candidates.first().cloned(),
            candidates
        }
    }
}