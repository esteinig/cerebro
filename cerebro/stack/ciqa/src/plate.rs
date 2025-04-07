
use std::{collections::{HashMap, HashSet}, path::Path};

use cerebro_gp::gpt::ClinicalContext;
use plotters::{coord::Shift, prelude::*};
use serde::{Deserialize, Serialize};
use colored::{ColoredString, Colorize};
use crate::{error::CiqaError, utils::read_tsv};

pub trait FromSampleType {
    fn from_sample_type(sample_type: SampleType) -> ClinicalContext;
}

impl FromSampleType for ClinicalContext {
    fn from_sample_type(sample_type: SampleType) -> ClinicalContext {
        match sample_type {
            SampleType::Csf => ClinicalContext::Csf,
            SampleType::Eye => ClinicalContext::Eye,
            _ => unimplemented!("Variant of sample type not implemented for FromSampleType -> Clinical Context")
        }
    }
}


// Define the PlateWell struct.
pub struct PlateWell {
    pub id: String,
    pub color: String,
    pub size: f64,
}

pub fn plot_plate() -> Result<(), CiqaError> {

    let rows = 8;
    let cols = 12;
    let margin = 50;
    let canvas_size = (800, 600);
    
    // Create PNG output using BitMapBackend.
    let png_backend = BitMapBackend::new("well_plate.png", canvas_size);
    let png_area = png_backend.into_drawing_area();
    draw_plate(png_area, rows, cols, margin)?;

    // Create SVG output using SVGBackend.
    let svg_backend = SVGBackend::new("well_plate.svg", canvas_size);
    let svg_area = svg_backend.into_drawing_area();
    draw_plate(svg_area, rows, cols, margin)
    
}
    

fn draw_plate<DB: DrawingBackend>(
    root: DrawingArea<DB, Shift>,
    rows: u32,
    cols: u32,
    margin: u32
) -> Result<(), CiqaError> where DB::ErrorType: 'static {

    // Fill the background with white.
    root.fill(&WHITE)?;

    // Get drawing area dimensions.
    let (width, height) = root.dim_in_pixel();
    let drawing_width = width - 2 * margin;
    let drawing_height = height - 2 * margin;
    let cell_width = drawing_width / cols;
    let cell_height = drawing_height / rows;

    // Create a vector of PlateWell instances.
    let mut wells: Vec<PlateWell> = Vec::new();
    let default_size = 1.0; // scaling factor
    let default_color = "blue".to_string();
    for row in 0..rows {
        for col in 0..cols {
            // Create unique ID e.g., "A1", "A2", … "H12"
            let id = format!("{}{}", (b'A' + row as u8) as char, col + 1);
            wells.push(PlateWell {
                id,
                color: default_color.clone(),
                size: default_size,
            });
        }
    }

    // Draw each well as a circle.
    for (i, well) in wells.iter().enumerate() {
        let row = i / cols as usize;
        let col = i % cols as usize;
        // Compute the center coordinates of the cell.
        let x_center = margin + col as u32 * cell_width + cell_width / 2;
        let y_center = margin + row as u32 * cell_height + cell_height / 2;

        // Determine the circle radius.
        let base_radius = (cell_width.min(cell_height) as f64) * 0.3;
        let radius = (base_radius * well.size).round() as u32;

        // Map the well's color string to a Plotters color.
        let plot_color = match well.color.as_str() {
            "blue" => BLUE,
            "red" => RED,
            "green" => GREEN,
            _ => BLACK,
        };

        // Draw the well as a filled circle.
        root.draw(&Circle::new(
            (x_center as i32, y_center as i32),
            radius,
            ShapeStyle::from(&plot_color).filled(),
        ))?;
    }

    // Finalize the drawing.
    root.present()?;

    Ok(())

}


// This enum represents the possible test results.
#[derive(Serialize, Deserialize, Debug, PartialEq, Clone)]
#[serde(rename_all = "lowercase")]
pub enum TestResult {
    Positive,
    Negative,
}

// This enum represents the possible test results.
#[derive(Serialize, Deserialize, Debug, PartialEq, Clone, clap::ValueEnum)]
#[serde(rename_all = "UPPERCASE")]
pub enum SampleType {
    Csf,
    Eye,
    Ntc,
    Pos,
    Env
}

// Represents a single test within an orthogonal method.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Test {
    pub test: Vec<String>,
    pub result: TestResult,
    pub taxa: Vec<String>,
    pub note: Option<String>,
}

// Represents an orthogonal method (for example, "culture", "pcr", etc.) along with its tests.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Orthogonal {
    pub method: String,
    pub detail: String,
    pub tests: Vec<Test>,
}

// Represents a single sample from the plate.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SampleReference {
    pub sample_id: String,
    pub sample_type: SampleType,
    pub result: Option<TestResult>,
    pub note: Option<String>,
    pub orthogonal: Vec<Orthogonal>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SampleReview {
    pub sample_id: String,
    pub result: Option<TestResult>,
    pub pathogen: Option<String>,
    pub note: Option<String>
}

// How to handle missing orthogonal tests for a reference sample when the 
// review call is positive
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, clap::ValueEnum)]
pub enum MissingOrthogonal {
    Indeterminate,
    ResultOnly        // not implemented yet, evaluate a positive against the overall reference result, not orthogonal testing
}

// --- Diagnostic outcome types ---

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, clap::ValueEnum)]
pub enum DiagnosticOutcome {
    TruePositive,
    FalsePositive,
    TrueNegative,
    FalseNegative,
    Indeterminate,
    NotConsidered,
}
impl std::fmt::Display for DiagnosticOutcome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            DiagnosticOutcome::TruePositive => "True Positive",
            DiagnosticOutcome::FalsePositive => "False Positive",
            DiagnosticOutcome::TrueNegative => "True Negative",
            DiagnosticOutcome::FalseNegative => "False Negative",
            DiagnosticOutcome::Indeterminate => "Indeterminate",
            DiagnosticOutcome::NotConsidered => "Not Considered",
        };
        write!(f, "{}", s)
    }
}

impl DiagnosticOutcome {
    /// Returns the string representation of the outcome with colored formatting.
    pub fn colored(&self) -> ColoredString {
        match self {
            Self::TruePositive | Self::TrueNegative => format!("{}", self).green(),
            Self::FalsePositive | Self::FalseNegative => format!("{}", self).red(),
            Self::Indeterminate => format!("{}", self).yellow(),
            Self::NotConsidered => format!("{}", self).white()
        }
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct DiagnosticReview {
    pub sample_id: String,
    pub outcome: DiagnosticOutcome,
    pub reference: SampleReference,
    pub review: Option<SampleReview>
}

/// Structure to hold computed diagnostic statistics.
#[derive(Serialize, Deserialize, Debug)]
pub struct DiagnosticStats {
    pub sensitivity: f64,
    pub specificity: f64,
    pub ppv: f64,
    pub npv: f64,
    pub total: usize
}

// The PlateReference bundles the entire list of samples.
#[derive(Serialize, Deserialize, Debug)]
pub struct ReferencePlate {
    pub reference: Vec<SampleReference>,
    pub review: Option<Vec<SampleReview>>,
    pub negative_controls: Vec<String>,
    pub samples: Vec<String>,
    pub missing_orthogonal: MissingOrthogonal
}

impl ReferencePlate {

    /// Create a new PlateReference by reading the reference from a JSON file and the optional review from a TSV file.
    pub fn new(reference_path: &Path, review_path: Option<&Path>, missing_orthogonal: MissingOrthogonal) -> Result<Self, CiqaError> {

        // Read the JSON file to get the vector of SampleReference.
        let reference_data = std::fs::read_to_string(reference_path)?;
        let reference: Vec<SampleReference> = serde_json::from_str(&reference_data)?;
        let negative_controls = Self::get_negative_controls(&reference);
        let samples = Self::get_samples(&reference);

        Ok(ReferencePlate { 
            reference, 
            review: if let Some(path) = review_path { Some(read_tsv(path, false, true)?) } else { None }, 
            negative_controls,
            samples,
            missing_orthogonal
        })
    }
    pub fn get_sample_reference(&self, sample_id: &str) -> Option<SampleReference> {
        self.reference
            .iter()
            .find(|r| r.sample_id == sample_id)
            .cloned()
    }
    pub fn get_negative_controls(reference: &Vec<SampleReference>) -> Vec<String> {
        reference
            .into_iter()
            .filter(|r| r.sample_type == SampleType::Ntc || r.sample_type == SampleType::Env)
            .map(|r| r.sample_id.clone())
            .collect()
    }
    pub fn get_samples(reference: &Vec<SampleReference>) -> Vec<String> {
        reference
            .into_iter()
            .filter(|r| r.sample_type == SampleType::Eye || r.sample_type == SampleType::Csf)
            .map(|r| r.sample_id.clone())
            .collect()
    }
    pub fn set_none(&mut self, sample_id: &Vec<String>) -> Result<(), CiqaError> {
        let ids: HashSet<&String> = sample_id.into_iter().collect();
        for sample in &mut self.reference {
            if ids.contains(&sample.sample_id) {
                sample.result = None;
            }
        }
        Ok(())
    }
    pub fn subset_sample_type(&mut self, sample_type: SampleType) -> Result<(), CiqaError> {
        self.reference = self.reference.iter().filter(|r| r.sample_type == sample_type).cloned().collect();
        Ok(())
    }
    /// Compute diagnostic outcomes by comparing each SampleReference with its matching SampleReview
    pub fn compute_diagnostic_review(&self) -> Result<Vec<DiagnosticReview>, CiqaError> {
        if let Some(ref reviews) = self.review {

            // Create a lookup by sample_id.
            let review_map: HashMap<&str, &SampleReview> = reviews
                .iter()
                .map(|r| (r.sample_id.as_str(), r))
                .collect();

            Ok(self.reference.iter().map(|reference| {
                
                if reference.sample_id == "DW-63-V64" {
                    log::info!("DW-63-V64")
                }

                let (outcome, review) = if let Some(review) = review_map.get(reference.sample_id.as_str()) {
                    (compare_sample_review(reference, review, self.missing_orthogonal.clone()), Some(*review))
                } else {
                    (DiagnosticOutcome::Indeterminate, None)
                };
                DiagnosticReview {
                    sample_id: reference.sample_id.clone(),
                    outcome,
                    reference: reference.clone(),
                    review: review.cloned()
                }
            }).collect())

        } else {
            log::warn!("No review data table was provided for diagnostic evaluation");
            Ok(Vec::new())
        }
    }


    /// Given a vector of DiagnosticReview, compute sensitivity, specificity,
    /// negative predictive value and positive predictive value, excluding any
    /// DiagnosticOutcome::Indeterminate results.
    pub fn compute_statistics(&self, diagnostic_review: Vec<DiagnosticReview>) -> DiagnosticStats {

        // Filter out DiagnosticOutcome::Indeterminate
        let filtered: Vec<DiagnosticReview> = diagnostic_review.into_iter()
            .filter(|d| !(d.outcome == DiagnosticOutcome::Indeterminate || d.outcome == DiagnosticOutcome::NotConsidered))
            .collect();

        let tp = filtered.iter().filter(|d| matches!(d.outcome, DiagnosticOutcome::TruePositive)).count();
        let fp = filtered.iter().filter(|d| matches!(d.outcome, DiagnosticOutcome::FalsePositive)).count();
        let tn = filtered.iter().filter(|d| matches!(d.outcome, DiagnosticOutcome::TrueNegative)).count();
        let fn_ = filtered.iter().filter(|d| matches!(d.outcome, DiagnosticOutcome::FalseNegative)).count();

        let sensitivity = if tp + fn_ > 0 {
            tp as f64 / (tp + fn_) as f64
        } else {
            0.0
        };

        let specificity = if tn + fp > 0 {
            tn as f64 / (tn + fp) as f64
        } else {
            0.0
        };

        let ppv = if tp + fp > 0 {
            tp as f64 / (tp + fp) as f64
        } else {
            0.0
        };

        let npv = if tn + fn_ > 0 {
            tn as f64 / (tn + fn_) as f64
        } else {
            0.0
        };

        let total: usize = filtered.len();

        DiagnosticStats {
            sensitivity,
            specificity,
            ppv,
            npv,
            total
        }
    }
    
    /// Reads JSON data from the given file path and deserializes it into a PlateReference.
    pub fn from_json<P: AsRef<Path>>(path: P) -> Result<Self, CiqaError> {
        let data = std::fs::read_to_string(path)?;
        let plate = serde_json::from_str::<ReferencePlate>(&data)?;
        Ok(plate)
    }

    /// Serializes the PlateReference into pretty-printed JSON and writes it to the given file path.
    pub fn to_json<P: AsRef<Path>>(&self, path: P) -> Result<(), CiqaError> {
        let data = serde_json::to_string_pretty(&self)?;
        std::fs::write(path, data)?;
        Ok(())
    }
}



/// Compare a single SampleReference and SampleReview and return a diagnostic outcome.
fn compare_sample_review(reference: &SampleReference, review: &SampleReview, missing_orthogonal: MissingOrthogonal) -> DiagnosticOutcome {

    // Independent of review result, if the reference result is missing assign indeterminate outcome
    if let None = reference.result {
        return DiagnosticOutcome::NotConsidered
    }

    match review.result {
        Some(TestResult::Positive) => {

            // For a positive review, we use the pathogen for comparison.
            if let Some(ref pathogen) = review.pathogen {

                // Positive review but no orthogonal tests conducted against which to assess <====== how should we handle these cases?
                if reference.orthogonal.is_empty() {
                    return match missing_orthogonal {
                        MissingOrthogonal::Indeterminate => DiagnosticOutcome::Indeterminate,
                        MissingOrthogonal::ResultOnly => match reference.result {
                            Some(TestResult::Negative) => DiagnosticOutcome::FalsePositive,
                            Some(TestResult::Positive) => DiagnosticOutcome::TruePositive,  // in the absence of orthogonal tests any pathogen in review will be a true positive
                            None => DiagnosticOutcome::Indeterminate,
                        }
                    }
                }
                
                // Look for any positive orthogonal test in the SampleReference that has a matching test string.
                for ortho in &reference.orthogonal {
                    for test in &ortho.tests {
                        if test.result == TestResult::Positive {
                            if test.taxa.iter().any(|t| t == pathogen) {
                                return DiagnosticOutcome::TruePositive;
                            }
                        }
                    }
                }
                // If no positive test in the reference matches the pathogen, then it's a false positive.
                DiagnosticOutcome::FalsePositive
            } else {
                DiagnosticOutcome::Indeterminate
            }
        }
        Some(TestResult::Negative) => {

            // When review is negative, simply compare the overall sample result.
            match reference.result {
                Some(TestResult::Negative) => DiagnosticOutcome::TrueNegative,
                Some(TestResult::Positive) => DiagnosticOutcome::FalseNegative,
                None => DiagnosticOutcome::Indeterminate,
            }
        }
        None => DiagnosticOutcome::Indeterminate,
    }
}