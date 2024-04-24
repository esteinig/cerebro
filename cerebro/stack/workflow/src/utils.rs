use colored::Colorize;
use std::{ffi::OsStr, path::{Path, PathBuf}};
use tabled::{Table, Alignment, Modify, Panel, object::Rows};
use needletail::parser::SequenceRecord;
use std::str::from_utf8;
use env_logger::Builder;
use env_logger::fmt::Color;
use log::{LevelFilter, Level};
use std::io::Write;

use crate::{
    sample::WorkflowSample, 
    quality::QualityControlSummary, 
    error::WorkflowError, 
    taxon::TaxonomyWarning
};

// Utility function to extract the database name as valid UTF-8
pub fn get_file_stem(file: &PathBuf) -> Result<String, WorkflowError> {
    match file.file_stem() {
        Some(name) => Ok(name
            .to_os_string()
            .into_string()
            .map_err(|_| WorkflowError::InvalidReferencePath)?),
        None => {
            return Err(WorkflowError::ReferenceNameExtraction(format!(
                "{:?}",
                file
            )))
        }
    }
}


// Convenience function to generate a QC table from multiple sample JSON
pub fn create_qc_table(samples: Vec<PathBuf>, output: &PathBuf, header: bool, ercc_input_mass: Option<f64>) -> Result<(), WorkflowError> {

    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(header).from_path(&output).unwrap();

    for file in samples {
        let sample = WorkflowSample::read_json(&file).expect(&format!("Failed to parse sample file: {}", file.display()));
        let row = QualityControlSummary::from(&sample.qc_module, None, ercc_input_mass, None)?;
        writer.serialize(&row).map_err(|err|WorkflowError::SerializeTableRow(err))?;     
    }
    Ok(())
}

/// Quality control summary table for WorkflowSample and WorkflowSampleSet
pub fn _qc_table(rows: Vec<QualityControlSummary>, file: Option<PathBuf>, header: bool) -> Result<(), WorkflowError> {

    match file {
        Some(file) => { 
            let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(header).from_path(&file).unwrap();

            for row in rows {
                writer.serialize(&row).unwrap();
            }
        },
        None => {
            // let mut table = Table::new(rows);
            // let header_color = tabled::color::Color::try_from(" ".blue().to_string()).unwrap();
            // let header_text = ansi_term::Style::new().bold().paint(format!("Quality Control")).to_string();
            // table.with(Panel::header(header_text))
            //     .with(Modify::new(Rows::first()).with(Alignment::center()).with(header_color));
            // println!("{}", table);
        }
    }
    Ok(())
}

pub fn _taxonomy_warning_table(warnings: Vec<&TaxonomyWarning>) -> Result<(), WorkflowError> {

    let mut table = Table::new(warnings);
    let header_color = tabled::color::Color::try_from(" ".red().to_string()).unwrap();
    let header_text = ansi_term::Style::new().bold().paint(format!("Reference Taxonomy Warnings")).to_string();
    let footer_text = ansi_term::Style::new().bold().paint(String::from(
        "These warnings are caused by mismatches between database taxonomic identifiers and the reference taxonomy. Records are not included in results."
    )).to_string();
    table.with(Panel::header(header_text))
         .with(Panel::footer(footer_text))
         .with(Modify::new(Rows::first()).with(Alignment::center()).with(header_color));

        println!("{}", table);

    Ok(())
}

pub fn get_colored_string(value: &str, color: &str) -> String {
    match color.to_lowercase().as_str() {
        "blue" => value.blue().to_string(),
        "green"  => value.green().to_string(),
        "yellow"  => value.yellow().to_string(),
        "red" => value.red().to_string(),
        "cyan" => value.cyan().to_string(),
        "magenta" => value.magenta().to_string(),
        _ => value.white().to_string()
    }
}

pub fn get_seq_record_identifier(rec: &SequenceRecord) -> Result<String, WorkflowError> {
    match from_utf8(rec.id())?.split(' ').next() {
        Some(rec_id) => Ok(rec_id.to_string()),
        None => return Err(WorkflowError::RecordIdentifierNotParsed)
    }
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

pub fn get_compression_writer(
    output: &std::path::PathBuf,
    output_format: &Option<niffler::compression::Format>, 
    compression_level: &Option<niffler::compression::Level>
) -> Result<Box<dyn std::io::Write>, WorkflowError> {

    let file_handle = std::io::BufWriter::new(std::fs::File::create(&output)?);

    let fmt = match output_format {
        None => niffler::Format::from_path(&output),
        Some(format) => *format,
    };

    let level = match compression_level {
        Some(level) => *level,
        None => niffler::compression::Level::Six
    };
    
    Ok(niffler::get_writer(Box::new(file_handle), fmt, level)?)

}


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

