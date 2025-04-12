
use std::{collections::{HashMap, HashSet}, fs::File, io::{BufRead, BufReader, BufWriter}, ops::Deref, path::{Path, PathBuf}};

use cerebro_gp::gpt::{ClinicalContext, Diagnosis, DiagnosticResult};
use plotters::{coord::Shift, prelude::*, style::{full_palette::GREY, text_anchor::{HPos, Pos, VPos}}};
use serde::{Deserialize, Serialize};
use colored::{ColoredString, Colorize};
use crate::{error::CiqaError, terminal::ReviewArgs, utils::{get_file_component, read_tsv, FileComponent}};
use rand::Rng;

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
    
    let png_backend = BitMapBackend::new("well_plate.png", canvas_size);
    let png_area = png_backend.into_drawing_area();
    draw_plate(png_area, rows, cols, margin)?;

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

    root.fill(&WHITE)?;

    let (width, height) = root.dim_in_pixel();

    let drawing_width = width - 2 * margin;
    let drawing_height = height - 2 * margin;

    let cell_width = drawing_width / cols;
    let cell_height = drawing_height / rows;

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
        let y_center = margin as i32 + row as i32 * cell_height as i32 + (cell_height / 2) as i32;

        let base_radius = (cell_width.min(cell_height) as f64) * 0.3;
        let radius = (base_radius * well.size).round() as u32;

        let plot_color = match well.color.as_str() {
            "blue" => BLUE,
            "red" => RED,
            "green" => GREEN,
            _ => BLACK,
        };

        root.draw(&Circle::new(
            (x_center as i32, y_center),
            radius,
            ShapeStyle::from(&plot_color).filled(),
        ))?;
    }

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
impl SampleReview {
    pub fn from_gpt(path: &Path) -> Result<Self, CiqaError> {
        
        let diagnostic_result = DiagnosticResult::from_json(path)?;

        let stem = path
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or(CiqaError::FileStemNotFound(path.to_path_buf()))?;

        let mut stem_parts = stem.split('.');

        let sample_id = stem_parts.next().ok_or(CiqaError::SampleIdentifierNotFound)?;
        let gpt_model = stem_parts.next().ok_or(CiqaError::ModelNameNotFound)?;
        
        let test_result =  if diagnostic_result.diagnosis == Diagnosis::Infectious || diagnostic_result.diagnosis == Diagnosis::InfectiousReview {
            Some(TestResult::Positive)
        } else if diagnostic_result.diagnosis == Diagnosis::NonInfectious || diagnostic_result.diagnosis == Diagnosis::NonInfectiousReview {
            Some(TestResult::Negative)
        } else {
            log::warn!("Diagnosis from MetaGPT is: {:?} ({})", diagnostic_result, sample_id);
            None // Diagnosis::Tumor
        };

        let pathogen_str = match diagnostic_result.pathogen {
            Some(pathogen_str) => {

                // First, perform the basic replacements and trimming.
                let pathogen = pathogen_str
                    .replace("'", "")
                    .replace("_", " ")
                    .trim_start_matches("s  ")
                    .to_string();

                let pathogen_parts: Vec<&str> = pathogen.split_whitespace().collect();

                // Check if there is at least one element to drop the last item (species variant in GTDB).
                let species = if pathogen_parts.len() > 2 {
                    pathogen_parts[..pathogen_parts.len()-1].join(" ")
                } else {
                    pathogen.to_string()
                };
    
                Some(format!("s__{}", species))
            },
            None => None
        };
        
        
        Ok(Self {
            sample_id: sample_id.to_string(),
            result: test_result,
            pathogen: pathogen_str,
            note: Some(format!("MetaGPT model: {gpt_model}"))
        })
    }
}


// Ensure CiqaError implements From<std::io::Error> or change the error handling accordingly.
pub fn read_all_sample_reviews(directory: &Path) -> Result<Vec<SampleReview>, CiqaError> {
    let mut sample_reviews = Vec::new();

    // Iterate over the entries in the provided directory
    for entry in std::fs::read_dir(directory)? {
        let entry = entry?;
        let path = entry.path();

        // Only process files that end with .json
        if path.is_file() {
            if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
                if ext.eq_ignore_ascii_case("json") {
                    // Use SampleReview::from_gpt to convert the JSON file into a SampleReview instance.
                    let review = SampleReview::from_gpt(&path)?;
                    sample_reviews.push(review);
                }
            }
        }
    }

    Ok(sample_reviews)
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
    Control,
    Unknown
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
            DiagnosticOutcome::Control => "Control",
            DiagnosticOutcome::Unknown => "Unknown",
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
            Self::Control => format!("{}", self).blue(),
            Self::NotConsidered => format!("{}", self).white(),
            Self::Unknown => format!("{}", self).purple(),

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
    pub name: String,
    pub sensitivity: f64,
    pub specificity: f64,
    pub ppv: f64,
    pub npv: f64,
    pub total: usize
}
impl DiagnosticStats {
    pub fn percent(&self) -> Self {
        Self {
            name: self.name.to_string(),
            sensitivity: self.sensitivity*100.0,
            specificity: self.specificity*100.0,
            ppv: self.ppv*100.0,
            npv: self.npv*100.0,
            total: self.total
        }
    }
}

/// Extension trait for `Vec<DiagnosticStats>` to write JSON to a file.
pub trait DiagnosticStatsVecExt {
    /// Writes this vector of diagnostic stats to the given path in JSON format.
    /// Returns an error if the file cannot be written or serialization fails.
    fn to_json(&self, path: &Path) -> Result<(), CiqaError>;
}

impl DiagnosticStatsVecExt for Vec<DiagnosticStats> {
    fn to_json(&self, path: &Path) -> Result<(), CiqaError> {
        let writer = BufWriter::new(
            File::create(path)?
        );
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }
}

pub trait RgbColor {
    fn from_hex_str(hex: &str) -> Result<Self, String> where Self: Sized;
}

impl RgbColor for RGBColor {
    fn from_hex_str(hex: &str) -> Result<Self, String> where Self: Sized {

        // Remove the leading '#' if it exists.
        let hex = hex.trim_start_matches('#');

        // Ensure the string is exactly 6 characters long
        if hex.len() != 6 {
            return Err(format!("Invalid hex code length for '{}'. Expected 6 characters.", hex));
        }

        // Parse each pair of hex digits into its corresponding u8 value.
        let r = u8::from_str_radix(&hex[0..2], 16)
            .map_err(|e| format!("Invalid red component in '{}': {}", hex, e))?;
        let g = u8::from_str_radix(&hex[2..4], 16)
            .map_err(|e| format!("Invalid green component in '{}': {}", hex, e))?;
        let b = u8::from_str_radix(&hex[4..6], 16)
            .map_err(|e| format!("Invalid blue component in '{}': {}", hex, e))?;

        Ok(RGBColor(r, g, b))
    }
}

pub struct Palette {
    pub colors: Vec<RGBColor>,
}

impl Palette {
    /// Constructs a Palette from a vector of hex color strings.
    /// Panics if any of the hex strings is invalid.
    pub fn from_hex(hex: Vec<String>) -> Self {
        let colors = hex
            .into_iter()
            .map(|hex_str| {
                // Note: using unwrap() here. In a real-world scenario, you might want to handle errors more gracefully.
                RGBColor::from_hex_str(&hex_str)
                    .unwrap_or_else(|err| panic!("Error parsing color '{}': {}", hex_str, err))
            })
            .collect();
        Palette { colors }
    }
    
    /// Returns a Palette with a set of default Monet palette colors.
    pub fn monet() -> Self {
        let monet_hex = vec![
            "#4e6d58".to_owned(),
            "#749e89".to_owned(),
            "#abccbe".to_owned(),
            "#e3cacf".to_owned(),
            "#c399a2".to_owned(),
            "#9f6e71".to_owned(),
            "#41507b".to_owned(),
            "#7d87b2".to_owned(),
            "#c2cae3".to_owned(),
        ];
        Palette::from_hex(monet_hex)
    }
    /// Returns a Palette with a set of default Manet palette colors.
    pub fn manet() -> Self {
        let manet_hex = vec![
            "#3b2319".to_string(),
            "#80521c".to_string(),
            "#d29c44".to_string(),
            "#ebc174".to_string(),
            "#ede2cc".to_string(),
            "#7ec5f4".to_string(),
            "#4585b7".to_string(),
            "#225e92".to_string(),
            "#183571".to_string(),
            "#43429b".to_string(),
            "#5e65be".to_string(),
        ];
        Palette::from_hex(manet_hex)
    }
    /// Returns a Palette with a set of default Hiroshige palette colors.
    pub fn hiroshige() -> Self {
        let hiroshige_hex = vec![
            "#e76254".to_string(),
            "#ef8a47".to_string(),
            "#f7aa58".to_string(),
            "#ffd06f".to_string(),
            "#ffe6b7".to_string(),
            "#aadce0".to_string(),
            "#72bcd5".to_string(),
            "#528fad".to_string(),
            "#376795".to_string(),
            "#1e466e".to_string(),
        ];
        Palette::from_hex(hiroshige_hex)
    }

    /// Returns a Palette with a set of default Cassatt2 palette colors.
    pub fn cassatt2() -> Self {
        let cassatt2_hex = vec![
            "#2d223c".to_string(),
            "#574571".to_string(),
            "#90719f".to_string(),
            "#b695bc".to_string(),
            "#dec5da".to_string(),
            "#c1d1aa".to_string(),
            "#7fa074".to_string(),
            "#466c4b".to_string(),
            "#2c4b27".to_string(),
            "#0e2810".to_string(),
        ];
        Palette::from_hex(cassatt2_hex)
    }
}

pub fn get_diagnostic_stats(args: &ReviewArgs, reference_plate: &mut ReferencePlate, review_name: &str) -> Result<DiagnosticStats, CiqaError> {

    if let Some(sample_id) = &args.set_none {
        reference_plate.set_none(sample_id)?;
    }

    if let Some(sample_type) = &args.sample_type {
        reference_plate.subset_sample_type(sample_type.clone())?;
    }

    let diagnostic_review = reference_plate.compute_diagnostic_review()?;

    for dr in &diagnostic_review {
        log::info!("{} => {}", dr.sample_id, dr.outcome.colored())
    }

    let stats = reference_plate.compute_statistics(
        &review_name, 
        diagnostic_review
    );

    Ok(stats)
}

// The PlateReference bundles the entire list of samples.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ReferencePlate {
    pub reference: Vec<SampleReference>,
    pub review: Option<Vec<SampleReview>>,
    pub negative_controls: Vec<String>,
    pub samples: Vec<String>,
    pub missing_orthogonal: MissingOrthogonal
}

impl ReferencePlate {

    /// Create a new PlateReference by reading the reference from a JSON file and the optional review from a TSV file.
    pub fn new(reference_path: &Path, review_path: Option<&Path>, missing_orthogonal: MissingOrthogonal, diagnostic_agent: bool) -> Result<Self, CiqaError> {

        // Read the JSON file to get the vector of SampleReference.
        let reference_data = std::fs::read_to_string(reference_path)?;
        let reference: Vec<SampleReference> = serde_json::from_str(&reference_data)?;

        let negative_controls = Self::get_negative_controls(&reference);
        let samples = Self::get_samples(&reference);

        let review = if let Some(path) = review_path { 
            if diagnostic_agent {
                Some(read_all_sample_reviews(path)?)
            } else {
                Some(read_tsv(path, false, true)?) 
            }
        } else { 
            None 
        };


        Ok(ReferencePlate { 
            reference, 
            review, 
            negative_controls,
            samples,
            missing_orthogonal
        })
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
    pub fn get_positive_controls(reference: &Vec<SampleReference>) -> Vec<String> {
        reference
            .into_iter()
            .filter(|r| r.sample_type == SampleType::Pos)
            .map(|r| r.sample_id.clone())
            .collect()
    }
    pub fn get_samples(reference: &Vec<SampleReference>) -> Vec<String> {
        reference
            .into_iter()
            .filter(|r| !([SampleType::Ntc, SampleType::Env, SampleType::Pos].contains(&r.sample_type)))
            .map(|r| r.sample_id.clone())
            .collect()
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
                
                let (outcome, review) = if let Some(review) = review_map.get(reference.sample_id.as_str()) {
                    (compare_sample_review(reference, review, self.missing_orthogonal.clone()), Some(*review))
                } else {

                    if [SampleType::Ntc, SampleType::Env, SampleType::Pos].contains(&reference.sample_type) {
                        (DiagnosticOutcome::Control, None)
                    } else {
                        (DiagnosticOutcome::Indeterminate, None)
                    }
                };
                DiagnosticReview {
                    sample_id: reference.sample_id.clone(),
                    outcome,
                    reference: reference.clone(),
                    review: review.cloned()
                }
            }).collect())

        } else {
            log::warn!("Review data was not provided to `ReferencePlate` for diagnostic evaluation");
            Ok(Vec::new())
        }
    }


    /// Given a vector of DiagnosticReview, compute sensitivity, specificity,
    /// negative predictive value and positive predictive value, excluding any
    /// DiagnosticOutcome::Indeterminate results.
    pub fn compute_statistics(&self, name: &str, diagnostic_review: Vec<DiagnosticReview>) -> DiagnosticStats {

        // Filter out DiagnosticOutcome::Indeterminate
        let filtered: Vec<DiagnosticReview> = diagnostic_review.into_iter()
            .filter(|d| !([
                    DiagnosticOutcome::Indeterminate, 
                    DiagnosticOutcome::NotConsidered, 
                    DiagnosticOutcome::Control, 
                    DiagnosticOutcome::Unknown
                ].contains(&d.outcome))
            )
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
            name: name.to_string(),
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


// This is our function to aggregate reviews from multiple ReferencePlate structs.
// It returns a new ReferencePlate where the `review` field contains a consensus of all SampleReviews.
pub fn aggregate_reference_plates(plates: Vec<ReferencePlate>) -> ReferencePlate {
    // Create a HashMap to group reviews by sample_id.
    let mut review_groups: HashMap<String, Vec<SampleReview>> = HashMap::new();
    for plate in &plates {
        if let Some(reviews) = &plate.review {
            for review in reviews {
                // Group plates reviews from multiple plates by sample_id
                review_groups
                    .entry(review.sample_id.clone())
                    .or_default()
                    .push(review.clone());
            }
        }
    }

    // Now compute a consensus SampleReview for each sample id.
    let consensus_reviews: Vec<SampleReview> = review_groups
        .into_iter()
        .map(|(sample_id, reviews)| {
            let mut positive_count = 0;
            let mut negative_count = 0;
            // Count pathogens from positive votes.
            let mut pathogen_counts: HashMap<String, usize> = HashMap::new();

            for rev in &reviews {
                match rev.result {
                    Some(TestResult::Positive) => {
                        positive_count += 1;
                        if let Some(ref pathogen) = rev.pathogen {
                            *pathogen_counts.entry(pathogen.clone()).or_insert(0) += 1;
                        }
                    }
                    Some(TestResult::Negative) => {
                        negative_count += 1;
                    }
                    // We ignore None/other cases for the vote count
                    _ => {}
                }
            }

            // When there is a tie, or if positives outnumber negatives,
            // we choose positive; otherwise negative.
            let consensus_result = if positive_count >= negative_count {
                Some(TestResult::Positive)
            } else {
                Some(TestResult::Negative)
            };

            // For positive cases, select the majority pathogen, if present.
            let consensus_pathogen = if consensus_result == Some(TestResult::Positive) {
                // If no pathogen was found, then leave it as None.
                pathogen_counts
                    .into_iter()
                    .max_by_key(|&(_, count)| count)
                    .map(|(pathogen, _)| pathogen)
            } else {
                None
            };

            // Construct a new SampleReview for this sample id. The note field here is optional;
            // you might decide to include more information (for example, vote counts).
            SampleReview {
                sample_id,
                result: consensus_result,
                pathogen: consensus_pathogen,
                note: Some("Consensus review".to_string()),
            }
        })
        .collect();

    // We assume that all the ReferencePlate structs have the same reference, negative_controls,
    // samples and missing_orthogonal. So we can use the first one to build the new plate.
    let first_plate = plates
        .first()
        .expect("At least one reference plate must be provided");

    ReferencePlate {
        reference: first_plate.reference.clone(),
        review: Some(consensus_reviews),
        negative_controls: first_plate.negative_controls.clone(),
        samples: first_plate.samples.clone(),
        missing_orthogonal: first_plate.missing_orthogonal.clone(),
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
                            } else {
                                // No exact match to pathogen species if the reference is genus dereference the pathogen species to pathogen genus for evaluation
                                for taxon in &test.taxa {
                                    if taxon.starts_with("g__") {
                                        let pathogen_genus = format!(
                                            "g__{}", 
                                            pathogen.replace("s__", "")
                                                .split_whitespace()
                                                .collect::<Vec<&str>>()
                                                .first()
                                                .unwrap_or(&"")
                                            );
                                        if taxon == &pathogen_genus {
                                            return DiagnosticOutcome::TruePositive
                                        }
                                    }
                                }
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
        None => DiagnosticOutcome::Unknown,
    }
}


/// An enum to specify which two measures to plot.
#[derive(Debug, Clone, Copy, PartialEq, clap::ValueEnum)]
pub enum StatsMode {
    SensSpec,
    PpvNpv,
}

/// Loads multiple JSON files (each file containing a `Vec<DiagnosticStats>`).
/// The key of the HashMap is taken from each file's stem (the filename without extension).
/// E.g. if a file is named "LungCancer.json", the key is "LungCancer".
pub fn load_diagnostic_stats_from_files(
    paths: Vec<PathBuf>,
) -> Result<HashMap<String, Vec<DiagnosticStats>>, CiqaError> {
    let mut data_map = HashMap::new();

    for path in paths {
        let file_stem = get_file_component(&path, FileComponent::FileStem)?;
        let reader = BufReader::new(File::open(&path)?);

        let stats_vec: Vec<DiagnosticStats> = serde_json::from_reader(reader)?;
        data_map.insert(file_stem, stats_vec);
    }

    Ok(data_map)
}



/// Plot a horizontal strip-plot
pub fn plot_stripplot<B: DrawingBackend>(
    backend: B,
    data: &HashMap<String, Vec<DiagnosticStats>>,
    mode: StatsMode,
    ref1: Option<f64>,
    ref2: Option<f64>,
    col1: Option<&RGBColor>,
    col2: Option<&RGBColor>
) -> Result<(), CiqaError> {

    // Convert the categories to a sorted list to have consistent ordering
    let mut categories: Vec<String> = data.keys().cloned().collect();
    categories.sort();

    let n_cats = categories.len();
    
    let root_area = backend.into_drawing_area();

    root_area.fill(&WHITE).map_err(|e| CiqaError::StripPlotError(e.to_string()))?;

    let mut chart = ChartBuilder::on(&root_area)
        .margin(25)
        .x_label_area_size(40)
        .y_label_area_size(120) // more space for category names
        .build_cartesian_2d(0f64..100f64, 0f64..n_cats as f64)
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;

    chart
        .configure_mesh()
        .disable_mesh()
        .set_tick_mark_size(LabelAreaPosition::Bottom, 12)
        .y_desc("")
        .y_labels(0)
        .disable_y_axis()
        .x_desc("")
        .x_label_formatter(&|x_val: &f64| format!("{:.0}%", x_val))
        .x_label_style(
            ("monospace", 14)
            .into_font()
            .into_text_style(&root_area)
            .pos(Pos::new(HPos::Center, VPos::Bottom))
        )
        .draw()
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;

    // We place each category at integer y-values: 0, 1, 2, ...
    // We'll label them near the middle of each integer (i + 0.5).
    for (idx, cat_name) in categories.iter().enumerate() {

        let y_center = idx as f64 + 0.5;

        chart.draw_series(std::iter::once(
            Text::new(
                cat_name.clone(),
                (0.0, y_center),
                ("monospace", 14)
                    .into_font()
                    .into_text_style(&root_area)
                    .pos(Pos::new(HPos::Right, VPos::Center)),
            ),
        ))
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
    }

    let point_radius: u32 = 4;

    let color_1 = match col1 {
        Some(color) => color.to_owned(),
        None => RGBColor(210, 210, 210)
    };

    let color_2 = match col2 {
        Some(color) => color.to_owned(),
        None => RGBColor(120, 120, 120)
    };

    let color_3 = *Palette::cassatt2()
        .colors
        .get(5)
        .unwrap_or(&RGBColor(70, 70, 70));

    if let Some(ref1) = ref1 {
        chart.draw_series(DashedLineSeries::new(
            vec![(ref1, 0.0), (ref1, n_cats as f64)],
            5, /* size = length of dash */
            10, /* spacing */
            ShapeStyle {
                color: color_1.into(),
                filled: false,
                stroke_width: 2,
            }
        ))
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
    }

    if let Some(ref2) = ref2 {
        
        chart.draw_series(DashedLineSeries::new(
            vec![(ref2, 0.0), (ref2, n_cats as f64)],
            5, /* size = length of dash */
            10, /* spacing */
            ShapeStyle {
                color: color_2.into(),
                filled: false,
                stroke_width: 2,
            }
        ))
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
    }

    // We'll use a small random jitter so multiple points don't overlap exactly
    let mut rng = rand::rng();

    let mut measure1_shapes = Vec::new();
    let mut measure2_shapes = Vec::new();

    let mut consensus_shapes = Vec::new();

    for (idx, cat_name) in categories.iter().enumerate() {
        let stats_vec = match data.get(cat_name) {
            Some(v) => v,
            None => continue,
        };

        for stat in stats_vec {


            let (x1, x2) = match mode {
                StatsMode::SensSpec => (stat.sensitivity, stat.specificity),
                StatsMode::PpvNpv   => (stat.ppv, stat.npv),
            };

            // y offsets for measure1 vs measure2
            let base_1 = idx as f64 + 0.3;
            let base_2 = idx as f64 + 0.7;

            // small random vertical jitter: ± 0.04
            let jitter1: f64 = (rng.random::<f64>() - 0.5) * 0.08;
            let jitter2: f64 = (rng.random::<f64>() - 0.5) * 0.08;

            let y1 = base_1 + jitter1;
            let y2 = base_2 + jitter2;

            if stat.name == "consensus" {
                 // measure1 is drawn as two concentric circles (fill + black stroke)
                 consensus_shapes.push(Circle::new(
                    (x1, y1),
                    point_radius,
                    ShapeStyle {
                        color: color_3.into(),
                        filled: true,
                        stroke_width: 1,
                    }
                ));
                consensus_shapes.push(Circle::new(
                    (x1, y1),
                    point_radius,
                    ShapeStyle {
                        color: BLACK.into(),
                        filled: false,
                        stroke_width: 1,
                    }
                ));

                // measure2 is also drawn with fill + black stroke
                consensus_shapes.push(Circle::new(
                    (x2, y2),
                    point_radius,
                    ShapeStyle {
                        color: color_3.into(),
                        filled: true,
                        stroke_width: 1,
                    }
                ));
                consensus_shapes.push(Circle::new(
                    (x2, y2),
                    point_radius,
                    ShapeStyle {
                        color: BLACK.into(),
                        filled: false,
                        stroke_width: 1,
                    }
                ));
            } else {
                // measure1 is drawn as two concentric circles (fill + black stroke)
                measure1_shapes.push(Circle::new(
                    (x1, y1),
                    point_radius,
                    ShapeStyle {
                        color: color_1.into(),
                        filled: true,
                        stroke_width: 1,
                    }
                ));
                measure1_shapes.push(Circle::new(
                    (x1, y1),
                    point_radius,
                    ShapeStyle {
                        color: BLACK.into(),
                        filled: false,
                        stroke_width: 1,
                    }
                ));

                // measure2 is also drawn with fill + black stroke
                measure2_shapes.push(Circle::new(
                    (x2, y2),
                    point_radius,
                    ShapeStyle {
                        color: color_2.into(),
                        filled: true,
                        stroke_width: 1,
                    }
                ));
                measure2_shapes.push(Circle::new(
                    (x2, y2),
                    point_radius,
                    ShapeStyle {
                        color: BLACK.into(),
                        filled: false,
                        stroke_width: 1,
                    }
                ));
            }
        }
    }

    chart
        .draw_series(measure2_shapes)
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?
        .label(if mode == StatsMode::SensSpec {
            "Specificity"
        } else {
            "NPV"
        })
        .legend(move |(x, y)| Circle::new((x + 15, y - 7), 5, color_2.filled()));

    chart
        .draw_series(measure1_shapes)
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?
        .label(if mode == StatsMode::SensSpec {
            "Sensitivity"
        } else {
            "PPV"
        })
        .legend(move |(x, y)| Circle::new((x + 15, y - 7), 5, color_1.filled()));
    
    
    chart
        .draw_series(consensus_shapes)
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?
        .label("Consensus")
        .legend(move |(x, y)| Circle::new((x + 15, y - 7), 5, color_3.filled()));


    // Finally configure and draw the legend box
    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::UpperMiddle)
        .border_style(&WHITE)
        .background_style(WHITE.filled())
        .draw()
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;

    // Present the drawing
    root_area
        .present()
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;

    Ok(())
}
