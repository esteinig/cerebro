use std::{ffi::OsStr, fs::File, io::{BufWriter, Write}, path::Path};
use cerebro_model::api::{cerebro::{model::Cerebro, schema::TestResult}, training::model::TrainingResultRecord};
use csv::{Writer, WriterBuilder};
use env_logger::Builder;
use env_logger::fmt::Color;
use log::{LevelFilter, Level};
use niffler::get_writer;
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


pub fn get_tsv_writer(file: &Path, header: bool) -> Result<Writer<Box<dyn Write>>, HttpClientError> {
    
    let buf_writer = BufWriter::new(File::create(&file)?);

    let writer = get_writer(
        Box::new(buf_writer), 
        niffler::Format::from_path(file), 
        niffler::compression::Level::Six
    )?;

    let csv_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(header)
        .from_writer(writer);

    Ok(csv_writer)
}

pub fn write_tsv<T: Serialize>(data: &Vec<T>, file: &Path, header: bool) -> Result<(), HttpClientError> {

    let mut writer = get_tsv_writer(file, header)?;

    for value in data {
        // Serialize each value in the vector into the writer
        writer.serialize(&value)?;
    }

    // Flush and complete writing
    writer.flush()?;
    
    Ok(())
}


pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

/// Attempts to infer the compression type from the file extension.
/// If the extension is not known, then Uncompressed is returned.
impl CompressionExt for niffler::compression::Format {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self {
        let path = Path::new(p);
        match path.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => Self::Gzip,
            Some(Some("bz") | Some("bz2")) => Self::Bzip,
            Some(Some("lzma")) => Self::Lzma,
            _ => Self::No,
        }
    }
}


#[derive(Serialize)]
struct SummaryRow {
    cerebro_id: String,
    cerebro_name: String,
    sample_id: String,
    sample_tags: String,
    sample_group: Option<String>,
    sample_type: Option<String>,
    sample_date: Option<String>,

    workflow_id: String,
    workflow_name: String,
    workflow_pipeline: String,
    workflow_version: String,

    run_id: String,
    run_date: String,
}

impl From<&Cerebro> for SummaryRow {
    fn from(c: &Cerebro) -> Self {

        SummaryRow {
            cerebro_id: c.id.to_string(),
            cerebro_name: c.name.clone(),

            sample_id: c.sample.id.to_string(),
            sample_tags: c.sample.tags.iter().map(String::from).collect::<Vec<_>>().join("-"),
            sample_group: c.sample.sample_group.clone(),
            sample_type: c.sample.sample_type.clone(),
            sample_date: c.sample.sample_date.clone(),

            workflow_id: c.workflow.id.to_string(),
            workflow_name: c.workflow.name.clone(),
            workflow_pipeline: c.workflow.pipeline.clone(),
            workflow_version: c.workflow.version.clone(),

            run_id: c.run.id.to_string(),
            run_date: c.run.date.clone(),
        }
    }
}

/// Build and write the TSV summary for a list of models.
pub fn write_models_summary(models: &[Cerebro], outfile: &Path) -> Result<(), HttpClientError> {
    let rows: Vec<SummaryRow> = models.iter().map(SummaryRow::from).collect();
    write_tsv(&rows, outfile, true)
}