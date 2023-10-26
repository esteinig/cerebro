use colored::Colorize;
use std::path::PathBuf;
use tabled::{Table, Alignment, Modify, Panel, object::Rows};
use super::{sample::WorkflowSample, quality::QualityControlSummary, error::WorkflowError, taxon::TaxonomyWarning};

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
        let row = QualityControlSummary::from(&sample.qc_module, None, ercc_input_mass, None, None, None).unwrap();
        writer.serialize(&row).unwrap();
            
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
