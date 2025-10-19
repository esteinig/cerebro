
use std::{collections::{HashMap, HashSet}, fs::File, io::{BufReader, BufWriter}, path::{Path, PathBuf}};
use cerebro_client::client::CerebroClient;
use cerebro_model::api::{cerebro::schema::{CerebroIdentifierSchema, MetaGpConfig, PrefetchData, PrevalenceContaminationConfig, SampleType, TestResult}};
use cerebro_pipeline::{modules::quality::{QcStatus, QualityControlSummary}};
use regex::Regex;
use statrs::distribution::{ContinuousCDF, StudentsT};
use meta_gpt::gpt::{SampleContext, Diagnosis, DiagnosticResult};
use plotters::{coord::Shift, prelude::*, style::text_anchor::{HPos, Pos, VPos}};
use serde::{Deserialize, Serialize};
use colored::{ColoredString, Colorize};
use crate::{error::CiqaError, terminal::ReviewArgs, utils::{get_file_component, read_tsv, write_tsv, FileComponent}};
use rand::Rng;

pub fn parse_dir_components(dir_name: &str) -> Option<(String,String,String,u32)> {
    // examples:
    // qwen3-14b-q8-0_clinical_tiered_tiered-threshold_default_3
    // qwen3-7b-q4_0_noclinical_xxx_42
    let re = Regex::new(
        r"^qwen3-(\d+b)-(q\d+(?:[_-]\d+)?)_[^-_]+_(clinical|noclinical).*?_(\d+)$"
    ).unwrap();
    let caps = re.captures(dir_name)?;
    let params = caps.get(1)?.as_str().to_string();
    let quant  = caps.get(2)?.as_str().to_string();
    let clinical = caps.get(3)?.as_str().to_string();
    let replicate: u32 = caps.get(4)?.as_str().parse().ok()?;
    Some((params, quant, clinical, replicate))
}

#[derive(Serialize)]
pub struct VramRow {
    pub params: String,
    pub quant: String,
    pub clinical: String,
    pub replicate: u32,
    pub gpu_index: u32,
    pub peak_vram: u64, // bytes
    pub source_dir: String,
}

#[derive(Serialize)]
pub struct SecondsRow {
    pub params: String,
    pub quant: String,
    pub clinical: String,
    pub replicate: u32,
    pub file: String,     // sample-level bench file name
    pub seconds: f32,
    pub source_dir: String,
}

pub fn write_tsv_vram(rows: &[VramRow], out: &Path) -> anyhow::Result<()> {
    let file = std::fs::File::create(out)?;
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_writer(file);
    w.write_record(["params","quant","clinical","replicate","gpu_index","peak_vram","source_dir"])?;
    for r in rows {
        w.write_record(&[
            &r.params, &r.quant, &r.clinical,
            &r.replicate.to_string(),
            &r.gpu_index.to_string(),
            &r.peak_vram.to_string(),
            &r.source_dir
        ])?;
    }
    w.flush()?;
    Ok(())
}

pub fn write_tsv_seconds(rows: &[SecondsRow], out: &Path) -> anyhow::Result<()> {
    let file = std::fs::File::create(out)?;
    let mut w = csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_writer(file);
    w.write_record(["params","quant","clinical","replicate","file","seconds","source_dir"])?;
    for r in rows {
        w.write_record(&[
            &r.params, &r.quant, &r.clinical,
            &r.replicate.to_string(),
            &r.file,
            &format!("{}", r.seconds),
            &r.source_dir
        ])?;
    }
    w.flush()?;
    Ok(())
}

pub trait FromSampleType {
    fn from_sample_type(sample_type: SampleType) -> SampleContext;
}

impl FromSampleType for SampleContext {
    fn from_sample_type(sample_type: SampleType) -> SampleContext {
        match sample_type {
            SampleType::Csf => SampleContext::Csf,
            SampleType::Eye => SampleContext::Eye,
            SampleType::Lod => SampleContext::Spike,
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
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum OrganismDomain {
    Bacteria,
    Eukaryota,
    Viruses
}

// Represents a single sample from the plate.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SampleReference {
    pub sample_id: String,
    pub sample_type: SampleType,
    pub result: Option<TestResult>,
    pub note: Option<String>,
    pub clinical: Option<String>,
    pub domain: Option<OrganismDomain>,
    pub orthogonal: Vec<Orthogonal>,
}

impl SampleReference {
    pub fn positive_taxa(&self) -> Option<Vec<String>> {
        let mut set = HashSet::new();

        for orth in &self.orthogonal {
            for test in &orth.tests {
                if let Some(note) = &test.note {
                    if note.contains("exclude_test") {
                        continue;
                    }
                }
                if test.result == TestResult::Positive {
                    for t in &test.taxa {
                        set.insert(t.clone());
                    }
                }
            }
        }

        if set.is_empty() {
            None
        } else {
            Some(set.into_iter().collect())
        }
    }
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

                    let pathogen_parts: Vec<&str> = pathogen_str.split_whitespace().collect();

                    // Remove genus or species variants where genus + species
                    // name is longer than two items 
                    let species = if pathogen_parts.len() > 2 {
                        pathogen_parts[..pathogen_parts.len()-1].join(" ")
                    } else {
                        pathogen_str.to_string()
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

impl DiagnosticOutcome {
    fn index(&self) -> usize {
        match self {
            Self::TruePositive   => 0,
            Self::FalsePositive  => 1,
            Self::TrueNegative   => 2,
            Self::FalseNegative  => 3,
            Self::Indeterminate  => 4,
            Self::NotConsidered  => 5,
            Self::Control        => 6,
            Self::Unknown        => 7,
        }
    }
}


#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct DiagnosticReview {
    pub sample_id: String,
    pub outcome: DiagnosticOutcome,
    pub reference: SampleReference,
    pub review: Option<SampleReview>
}

/// Structure to hold computed diagnostic statistics.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct DiagnosticStats {
    pub name: String,
    pub sensitivity: f64,
    pub specificity: f64,
    pub ppv: f64,
    pub npv: f64,
    pub total: usize,
    pub true_positive: usize,
    pub true_negative: usize,
    pub false_positive: usize,
    pub false_negative: usize,
    pub data: Vec<DiagnosticReview>,
    pub data_filtered: Vec<DiagnosticReview>
}
impl DiagnosticStats {
    pub fn percent(&self) -> Self {
        Self {
            name: self.name.to_string(),
            sensitivity: self.sensitivity*100.0,
            specificity: self.specificity*100.0,
            ppv: self.ppv*100.0,
            npv: self.npv*100.0,
            total: self.total,
            true_positive: self.true_positive,
            true_negative: self.true_negative,
            false_positive: self.false_positive,
            false_negative: self.false_negative,
            data: self.data.clone(),
            data_filtered: self.data_filtered.clone()
        }
    }
}

impl std::fmt::Display for DiagnosticStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Name: {}  Sensitivity: {:.2}%  Specificity: {:.2}%  PPV: {:.2}%  NPV: {:.2}%  Total: {}",
            self.name,
            self.sensitivity,
            self.specificity,
            self.ppv,
            self.npv,
            self.total
        )
    }
}



#[derive(Serialize, Deserialize, Debug)]
pub struct ConsensusDiagnosticSummary {
    pub sensitivity: f64,
    pub specificity: f64,
    pub ppv: f64,
    pub npv: f64,
}
impl Default for ConsensusDiagnosticSummary {
    fn default() -> Self {
        Self {
            sensitivity: 0.0,
            specificity: 0.0,
            ppv: 0.0,
            npv: 0.0,
        }
    }
}
impl ConsensusDiagnosticSummary {
    pub fn from_stats(data: Vec<DiagnosticStats>) -> Self {
        // Find the entry with name == "consensus"
        if let Some(consensus) = data.into_iter().find(|s| s.name == "consensus") {
            Self {
                sensitivity: consensus.sensitivity,
                specificity: consensus.specificity,
                ppv: consensus.ppv,
                npv: consensus.npv,
            }
        } else {
            // If not found, return default (all zeros)
            Self::default()
        }
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct DiagnosticSummary {
    pub mean_sensitivity: f64,
    pub mean_specificity: f64,
    pub ci_sensitivity: Vec<f64>,
    pub ci_specificity: Vec<f64>,
    pub mean_ppv: f64,
    pub mean_npv: f64,
    pub ci_ppv: Vec<f64>,
    pub ci_npv: Vec<f64>
}
impl Default for DiagnosticSummary {
    fn default() -> Self {
        Self {
            mean_sensitivity: 0.0,
            mean_specificity: 0.0,
            ci_sensitivity: vec![0.0, 0.0],
            ci_specificity: vec![0.0, 0.0],
            mean_ppv: 0.0,
            mean_npv: 0.0,
            ci_ppv: vec![0.0, 0.0],
            ci_npv: vec![0.0, 0.0],
        }
    }
}
impl DiagnosticSummary {
    pub fn from_stats(data: Vec<DiagnosticStats>) -> Self {

        // Filter out any stats with name "consensus"
        let filtered: Vec<DiagnosticStats> = data
            .into_iter()
            .filter(|s| s.name != "consensus")
            .collect();
        
        // If no items remain, return zeros
        if filtered.is_empty() {
            return DiagnosticSummary::default()
        }
    
        // Extract each metric into its own vector
        let sensitivities: Vec<f64> = filtered.iter().map(|s| s.sensitivity).collect();
        let specificities: Vec<f64> = filtered.iter().map(|s| s.specificity).collect();
        let ppvs: Vec<f64> = filtered.iter().map(|s| s.ppv).collect();
        let npvs: Vec<f64> = filtered.iter().map(|s| s.npv).collect();
    
        // Compute mean and 95% CI for each metric
        let (mean_sensitivity, ci_sensitivity) = compute_mean_and_ci(&sensitivities);
        let (mean_specificity, ci_specificity) = compute_mean_and_ci(&specificities);
        let (mean_ppv, ci_ppv) = compute_mean_and_ci(&ppvs);
        let (mean_npv, ci_npv) = compute_mean_and_ci(&npvs);
    
        DiagnosticSummary {
            mean_sensitivity,
            ci_sensitivity,
            mean_specificity,
            ci_specificity,
            mean_ppv,
            mean_npv,
            ci_ppv,
            ci_npv,
        }
    }
}


/// Helper function to compute the mean and 95% confidence interval (CI) for a slice of f64 values.
/// The CI is computed using the t-distribution.
fn compute_mean_and_ci(values: &[f64]) -> (f64, Vec<f64>) {
    let n = values.len();
    if n == 0 {
        return (0.0, vec![0.0, 0.0]);
    }
    
    let mean: f64 = values.iter().sum::<f64>() / (n as f64);
    
    // If there is only one element, return a degenerate CI (both limits equal the mean).
    if n == 1 {
        return (mean, vec![mean, mean]);
    }
    
    // Compute sample variance (using n-1 for an unbiased estimate).
    let variance: f64 = values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / ((n - 1) as f64);
    let std_dev = variance.sqrt();
    let std_err = std_dev / (n as f64).sqrt();
    
    // Initialize a Student's t-distribution with n-1 degrees of freedom.
    let t_dist = StudentsT::new(0.0, 1.0, n as f64 - 1.0)
        .expect("Failed to create Student's t-distribution");
    
    // For a 95% CI, get the 97.5th percentile.
    let t_val = t_dist.inverse_cdf(0.975);
    
    let margin = t_val * std_err;
    let lower = mean - margin;
    let upper = mean + margin;
    
    (mean, vec![lower, upper])
}

pub fn average_replicate_certainty(
    data: &[Vec<DiagnosticReview>],
    reference: Option<&[DiagnosticReview]>,
) -> f64 {
    let nrows = if let Some(r) = reference {
        r.len()
    } else if !data.is_empty() {
        data[0].len()
    } else {
        return 0.0;
    };
    if nrows == 0 { return 0.0; }

    let mut total_good = 0usize;
    let mut total_count = 0usize;

    for row in 0..nrows {
        for col in data {
            let rev = &col[row];
            match rev.outcome {
                DiagnosticOutcome::TruePositive | DiagnosticOutcome::TrueNegative => { total_good += 1; total_count += 1; },
                DiagnosticOutcome::FalsePositive
                | DiagnosticOutcome::FalseNegative | DiagnosticOutcome::Indeterminate => { total_count += 1; },
                DiagnosticOutcome::NotConsidered
                | DiagnosticOutcome::Control
                | DiagnosticOutcome::Unknown => {}
            }
        }
    }

    if total_count == 0 {
        0.0
    } else {
        ((total_good as f64) / (total_count as f64))*100.0
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct DiagnosticData {
    pub consensus: ConsensusDiagnosticSummary,
    pub summary: DiagnosticSummary,
    pub stats: Vec<DiagnosticStats>
}
impl DiagnosticData {
    pub fn from(data: Vec<DiagnosticStats>) -> Self {
        Self {
            stats: data.clone(),
            summary: DiagnosticSummary::from_stats(data.clone()),
            consensus: ConsensusDiagnosticSummary::from_stats(data)
        }
    }
    pub fn to_json(&self, path: &Path) -> Result<(), CiqaError> {
        let writer = BufWriter::new(
            File::create(path)?
        );
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }
    pub fn from_json<P: AsRef<Path>>(path: P) -> Result<Self, CiqaError> {
        let data: String = std::fs::read_to_string(path)?;
        let diagnostic_data = serde_json::from_str::<DiagnosticData>(&data)?;
        Ok(diagnostic_data)
    }
    pub fn get_consensus_stats(&self) -> Option<DiagnosticStats> {
        self.stats.iter().find(|s| s.name == "consensus").cloned()
    }
    pub fn get_consensus_reviews_filtered(&self) -> Option<Vec<DiagnosticReview>> {
        if let Some(stats) = self.get_consensus_stats() {
            Some(stats.data_filtered)
        } else {
            None
        }
    }
    pub fn get_reference_from_json<P: AsRef<Path>>(path: P) -> Result<Option<Vec<DiagnosticReview>>, CiqaError> {
        let reference_data = Self::from_json(path)?;

        let reference_review: Vec<Vec<DiagnosticReview>> = reference_data.stats
            .iter()
            .filter_map(|d| {
                if d.name == String::from("consensus") { None } else { Some(d.data.clone()) }
            }).collect();
        
        Ok(reference_review.first().cloned())
    }
    pub fn plot_summary(&self, output: &Path, title: Option<&str>, width: u32, height: u32, reference: Option<PathBuf>, header_text: Option<&str>, consensus_stats: Option<&DiagnosticStats>) -> Result<(), CiqaError> {

        let data_columns: Vec<Vec<DiagnosticReview>> = self.stats
            .iter()
            .filter_map(|d| {
                if d.name == String::from("consensus") { None } else { Some(d.data.clone()) }
            }).collect();
        
        let consensus_column = self.stats
            .iter()
            .find_map(|d| {
                if d.name == String::from("consensus") { Some(d.data.clone()) } else { None }
            });
        
        let reference_column = match reference {
            Some(ref_path) => Self::get_reference_from_json(&ref_path)?,
            None => None
        };
        

        plot_diagnostic_matrix(
            &data_columns, 
            reference_column.as_deref(), 
            consensus_column.as_deref(), 
            output, 
            &Palette::diagnostic_review2(), 
            CellShape::Circle, 
            PanelColumnHeader::Panel,
            title,
            header_text,
            consensus_stats
        )
    }
    pub fn split_replicates_consensus(&self) -> (Vec<DiagnosticStats>, Vec<DiagnosticStats>) {
        let mut replicates = Vec::new();
        let mut consensus = Vec::new();

        for item in &self.stats {
            if item.name.to_lowercase().contains("consensus") {
                consensus.push(item.clone());
            } else {
                replicates.push(item.clone());
            }
        }

        (replicates, consensus)
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

pub enum PaletteName {
    Monet,
    Manet,
    Morgenstern,
    Hiroshige,
    Cassatt2,
    Paquin,
    DiagnosticReview,
    Custom(String)
}

pub struct Palette {
    pub name: PaletteName,
    pub colors: Vec<RGBColor>,
}

impl Palette {
    /// Constructs a Palette from a vector of hex color strings.
    /// Panics if any of the hex strings is invalid.
    pub fn from_hex(name: PaletteName, hex: Vec<String>) -> Self {
        let colors = hex
            .into_iter()
            .map(|hex_str| {
                // Note: using unwrap() here. In a real-world scenario, you might want to handle errors more gracefully.
                RGBColor::from_hex_str(&hex_str)
                    .unwrap_or_else(|err| panic!("Error parsing color '{}': {}", hex_str, err))
            })
            .collect();
        Palette { colors, name }
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
        Palette::from_hex(PaletteName::Monet, monet_hex)
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
        Palette::from_hex(PaletteName::Manet, manet_hex)
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
        Palette::from_hex(PaletteName::Hiroshige, hiroshige_hex)
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
        Palette::from_hex(PaletteName::Cassatt2, cassatt2_hex)
    }
    /// Returns a Palette with a set of default Cassatt2 palette colors.
    pub fn morgenstern() -> Self {
        let morgenstern_hex = vec![
            "#98768e".to_string(),
            "#b08ba5".to_string(),
            "#c7a2b6".to_string(),
            "#dfbbc8".to_string(),
            "#ffc680".to_string(),
            "#ffb178".to_string(),
            "#db8872".to_string(),
            "#a56457".to_string(),
        ];
        Palette::from_hex(PaletteName::Morgenstern, morgenstern_hex)
    }
    /// Returns a Palette with the "Paquin" palette colors.
    pub fn paquin() -> Self {
        let paquin_hex = vec![
            "#831818".to_owned(),
            "#c62320".to_owned(),
            "#f05b43".to_owned(),
            "#f78462".to_owned(),
            "#feac81".to_owned(),
            "#f7dea3".to_owned(),
            "#ced1af".to_owned(),
            "#98ab76".to_owned(),
            "#748f46".to_owned(),
            "#47632a".to_owned(),
            "#275024".to_owned(),
        ];
        Palette::from_hex(PaletteName::Paquin, paquin_hex)
    }
    /// Returns a Palette with a set of default colors for the diagnostic review plot
    pub fn diagnostic_review() -> Self {

        //  TruePositive   => 0,
        //  FalsePositive  => 1,
        //  TrueNegative   => 2,
        //  FalseNegative  => 3,
        //  Indeterminate  => 4,
        //  NotConsidered  => 5,
        //  Control        => 6,
        //  Unknown        => 7,

        let review_hex = vec![
            "#748f46".to_string(), 
            "#f05b43".to_string(),
            "#98ab76".to_string(),
            "#f78462".to_string(),
            "#f7dea3".to_string(),
            "#d3d3d3".to_string(),
            "#aadce0".to_string(),
            "#d3d3d3".to_string(),
        ];
        Palette::from_hex(PaletteName::DiagnosticReview, review_hex)
    }
    /// Returns a Palette with a set of default colors for the diagnostic review plot
    pub fn diagnostic_review2() -> Self {

        //  TruePositive   => 0,
        //  FalsePositive  => 1,
        //  TrueNegative   => 2,
        //  FalseNegative  => 3,
        //  Indeterminate  => 4,
        //  NotConsidered  => 5,
        //  Control        => 6,
        //  Unknown        => 7,

        let review_hex = vec![
            "#768f84".to_string(), 
            "#f78462".to_string(),
            "#B7C6C0".to_string(),
            "#fbc1b0".to_string(),
            "#dfbf88".to_string(),
            "#d3d3d3".to_string(),
            "#89a7c0".to_string(),
            "#d3d3d3".to_string(),
        ];
        Palette::from_hex(PaletteName::DiagnosticReview, review_hex)
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
        if dr.outcome == DiagnosticOutcome::FalseNegative || dr.outcome == DiagnosticOutcome::FalsePositive {
            log::info!("{} => {:<14}{}", dr.sample_id, dr.outcome.colored(), match dr.review { 
                Some(ref review) => match &review.pathogen {
                    Some(taxstr) => format!(" => {}", taxstr), 
                    None => "".to_string()
                }
                None => "".to_string() 
            })
        }
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
    pub fn from_path(reference_path: &Path) -> Result<Self, CiqaError> {

        let reference_data = std::fs::read_to_string(reference_path)?;
        let reference: Vec<SampleReference> = serde_json::from_str(&reference_data)?;

        let negative_controls = Self::get_negative_controls(&reference);
        let samples = Self::get_samples(&reference);

        
        Ok(ReferencePlate { 
            reference, 
            review: None, 
            negative_controls,
            samples,
            missing_orthogonal: MissingOrthogonal::Indeterminate
        })
    }
    /// Create a new PlateReference by reading the reference from a JSON file and the optional review from file (TSV or META-GPT).
    pub fn from_review(reference_path: &Path, review_path: Option<&Path>, missing_orthogonal: MissingOrthogonal, diagnostic_agent: bool) -> Result<Self, CiqaError> {

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
    /// DiagnosticOutcome::Indeterminate | NotConsidere | Control | Unknown.
    pub fn compute_statistics(&self, name: &str, diagnostic_review: Vec<DiagnosticReview>) -> DiagnosticStats {

        // Filter out DiagnosticOutcome::Indeterminate
        let filtered: Vec<DiagnosticReview> = diagnostic_review.iter()
            .filter(|d| !([
                    DiagnosticOutcome::Indeterminate, 
                    DiagnosticOutcome::NotConsidered, 
                    DiagnosticOutcome::Control, 
                    DiagnosticOutcome::Unknown
                ].contains(&d.outcome))
            )
            .cloned()
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
            total, 
            true_positive: tp,
            false_positive: fp,
            true_negative: tn,
            false_negative: fn_,
            data: diagnostic_review,
            data_filtered: filtered
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

    pub fn prevalence_contamination(
        &self, 
        client: &CerebroClient,
        contam_config: &PrevalenceContaminationConfig
    ) -> Result<HashMap<String, HashSet<String>>, CiqaError> {

        log::info!(
            "Fetching prevalence contamination taxids from database '{}' and project '{}'", 
            client.db.clone().unwrap_or_default(), client.project.clone().unwrap_or_default()
        );

        let dna_taxids = client.get_prevalence_contamination(contam_config, vec![String::from("DNA")])?;
        let rna_taxids = client.get_prevalence_contamination(contam_config, vec![String::from("RNA")])?;
        
        Ok(HashMap::from_iter(vec![
            ("DNA".to_string(), dna_taxids),
            ("RNA".to_string(), rna_taxids)
        ]))
    }

    pub fn prefetch(
        &self, 
        client: &CerebroClient,
        output: &Path, 
        config: &MetaGpConfig, 
        prevalence_contamination: HashMap<String, HashSet<String>>
    ) -> Result<PrefetchData, CiqaError> {
        
        log::info!("Fetching primary threshold data for sample '{}'", config.sample);

        let (primary, primary_contamination) = client.get_taxa( 
            &CerebroIdentifierSchema::from_gp_config(config), 
            &config.filter_configs.primary, 
            prevalence_contamination.clone(),
            config.contamination.outliers.primary
        )?;

        log::info!("Fetching secondary threshold data for sample '{}'", config.sample);


        let (secondary, secondary_contamination) = client.get_taxa( 
            &CerebroIdentifierSchema::from_gp_config(config), 
            &config.filter_configs.secondary, 
            prevalence_contamination.clone(),
            config.contamination.outliers.secondary
        )?;

        log::info!("Fetching target threshold data for sample '{}'", config.sample);
        
        let (target, target_contamination) = client.get_taxa( 
            &CerebroIdentifierSchema::from_gp_config(config), 
            &config.filter_configs.target, 
            prevalence_contamination.clone(),
            config.contamination.outliers.target
        )?;

        let mut prefetch_data = PrefetchData::new(
            primary, 
            secondary, 
            target, 
            primary_contamination,
            secondary_contamination,
            target_contamination,
            &config.filter_configs.primary,
            &config.filter_configs.secondary,
            &config.filter_configs.target,
            config
        );

        prefetch_data.prune();
        prefetch_data.to_json(output)?;

        Ok(prefetch_data)
        
    }

    pub fn write_tsv(&self, output: &Path, species_rank: bool) -> Result<(), CiqaError> {
        let rows = reference_plate_to_tsv_rows(&self, species_rank);
        write_tsv(&rows, output, true)
    }
}


#[derive(Serialize)]
struct PlateTsvRow {
    /// sample identifier
    sample_id: String,
    /// sample type
    sample_type: SampleType,
    /// overall result from reference
    result: Option<TestResult>,
    /// semicolon-separated list of taxa seen in any positive orthogonal test
    taxa: String,
    /// domain of positive taxa if provided (single one only for now)
    domain: Option<OrganismDomain>,
    /// organisms any of the orthogonal tested for 
    tested: String
}

// Build the TSV rows from a ReferencePlate
fn reference_plate_to_tsv_rows(plate: &ReferencePlate, species_rank: bool) -> Vec<PlateTsvRow> {
    plate.reference.iter().map(|r| {
        let mut taxa_vec: Vec<String> = Vec::new();
        let mut taxa_ref: Vec<String> = Vec::new();
        for ortho in &r.orthogonal {
            'tests: for test in &ortho.tests {
                if let Some(note) = &test.note {
                    if note.contains("exclude_test") {
                        continue 'tests;
                    }
                }
                // Collect positives
                if test.result == TestResult::Positive {
                    'taxa: for t in &test.taxa {
                        if species_rank && !t.starts_with("s__") {  
                            continue 'taxa;
                        }
                        if !taxa_vec.iter().any(|x| x == t) {
                            taxa_vec.push(t.clone());
                        }
                    }
                }
                // Collect references
                for tested in &test.test {
                    if species_rank && !tested.starts_with("s__")  {
                        continue 'tests;
                    }
                    if !taxa_ref.iter().any(|x| x == tested) {
                        taxa_ref.push(tested.clone());
                    }
                }
            }
        }
        let taxa_joined = taxa_vec.join(";");
        let refs_joined = taxa_ref.join(";");

        PlateTsvRow {
            sample_id: r.sample_id.clone(),
            sample_type: r.sample_type.clone(),
            result: r.result.clone(),
            taxa: taxa_joined,
            domain: r.domain.clone(),
            tested: refs_joined
        }
    }).collect()
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

    // Independent of review result, if the reference result is missing assign not considered outcome
    if let None = reference.result {
        if [SampleType::Ntc, SampleType::Env, SampleType::Pos].contains(&reference.sample_type) {
            return DiagnosticOutcome::Control
        } else {
            return DiagnosticOutcome::NotConsidered
        }
    }

    // Independent of review result, if the reference note contains "no_data"
    if reference.note.clone().is_some_and(|x| x.contains("no_data")) {
        match missing_orthogonal {
            MissingOrthogonal::Indeterminate => return DiagnosticOutcome::Indeterminate,
            MissingOrthogonal::ResultOnly => {}  // continue in the result only mode
        }
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

                        // Exclude any tests that have the "exclude_test" note usually where confirmatory testing
                        // with another test has not found the reference pathogen.
                        if let Some(note) = &test.note {
                            if note.contains("exclude_test") {
                                continue;
                            }
                        }

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
                        } else {
                            continue; // if test result is negative skip this comparison
                        }
                    }
                }
                // If no positive test in the reference matches the pathogen, it's a false positive.
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
/// The key of the HashMap is taken from each file's stem (the filename without extension)..
pub fn load_diagnostic_stats_from_files(
    paths: Vec<PathBuf>,
) -> Result<HashMap<String, DiagnosticData>, CiqaError> {
    let mut data_map = HashMap::new();

    for path in paths {
        let file_stem = get_file_component(&path, FileComponent::FileStem)?;
        let reader = BufReader::new(File::open(&path)?);

        let data: DiagnosticData = serde_json::from_reader(reader)?;
        data_map.insert(file_stem, data);
    }

    Ok(data_map)
}



/// Plot a horizontal strip-plot
pub fn plot_stripplot<B: DrawingBackend>(
    backend: B,
    data: &HashMap<String, DiagnosticData>,
    mode: StatsMode,
    ref1: Option<f64>,
    ref2: Option<f64>,
    col1: Option<&RGBColor>,
    col2: Option<&RGBColor>,
    col3: Option<&RGBColor>,
    ci: bool,
    file_order: Vec<PathBuf>,
    legend_position: Option<SeriesLabelPosition>,
    boxplot: bool,
    barplot: bool,
    y_labels: Option<Vec<String>>,  // <-- new argument
) -> Result<(), CiqaError> {

    let mut file_labels = Vec::new();
    for file in file_order {
        file_labels.push(
            get_file_component(
                &file, 
                FileComponent::FileStem
            )?
        )
    }

    // Convert the categories to a sorted list to have consistent ordering
    let mut categories: Vec<String> = data.keys().cloned().collect();

    if file_labels.is_empty() {
        categories.sort();
    } else {
        categories = file_labels;
    }; 

    // Decide which labels to use
    let labels: Vec<String> = match &y_labels {
        Some(custom) if custom.len() == categories.len() => custom.clone(),
        _ => categories.clone(), // fallback to sample names
    };

    let n_cats = labels.len();
    
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
    for (idx, label) in labels.iter().enumerate() {

        let y_center = idx as f64 + 0.5;

        chart.draw_series(std::iter::once(
            Text::new(
                label.clone(),
                (0.0, y_center),
                ("monospace", 14)
                    .into_font()
                    .into_text_style(&root_area)
                    .pos(Pos::new(HPos::Right, VPos::Center)),
            ),
        ))
        .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
    }


    // Helper: simple quantile interpolation
    let quantile = |sorted: &Vec<f64>, p: f64| -> f64 {
        let n = sorted.len();
        if n == 0 { return 0.0; }
        let idx = p * (n - 1) as f64;
        let lo = idx.floor() as usize;
        let hi = idx.ceil() as usize;
        if lo == hi {
            sorted[lo]
        } else {
            sorted[lo] + (sorted[hi] - sorted[lo]) * (idx - lo as f64)
        }
    };

    // Optionally draw boxplot overlays
    if boxplot {
        let box_style = ShapeStyle {
            color: RGBAColor(100, 100, 100, 0.2),
            stroke_width: 0,
            filled: true,
        };
        let whisker_style = ShapeStyle {
            color: BLACK.into(),
            stroke_width: 1,
            filled: false,
        };
        for (idx, cat_name) in categories.iter().enumerate() {
            if let Some(diag_data) = data.get(cat_name) {
                let stats = &diag_data.stats;
                // Measure 1 values
                let mut vals1: Vec<f64> = stats.iter().map(|s| match mode {
                    StatsMode::SensSpec => s.sensitivity,
                    StatsMode::PpvNpv   => s.ppv,
                }).collect();
                vals1.sort_by(|a, b| a.partial_cmp(b).unwrap());
                if !vals1.is_empty() {
                    let q1 = quantile(&vals1, 0.25);
                    let med = quantile(&vals1, 0.5);
                    let q3 = quantile(&vals1, 0.75);
                    let iqr = q3 - q1;
                    let lw = *vals1.iter()
                        .filter(|v| **v >= q1 - 1.5 * iqr)
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap();
                    let uw = *vals1.iter()
                        .filter(|v| **v <= q3 + 1.5 * iqr)
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap();
                    let center_y1 = idx as f64 + 0.3;
                    let y0 = center_y1 - 0.1;
                    let y1 = center_y1 + 0.1;
                    // IQR box
                    chart.draw_series(std::iter::once(
                        Rectangle::new([(q1, y0), (q3, y1)], box_style.clone())
                    ))
                    .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
                    // Median line
                    chart.draw_series(std::iter::once(
                        PathElement::new(vec![(med, y0), (med, y1)], whisker_style.clone())
                    ))
                    .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
                    // Whiskers
                    chart.draw_series(std::iter::once(
                        PathElement::new(vec![(lw, center_y1), (q1, center_y1)], whisker_style.clone())
                    ))
                    .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
                    chart.draw_series(std::iter::once(
                        PathElement::new(vec![(q3, center_y1), (uw, center_y1)], whisker_style.clone())
                    ))
                    .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
                }
                // Measure 2 values
                let mut vals2: Vec<f64> = stats.iter().map(|s| match mode {
                    StatsMode::SensSpec => s.specificity,
                    StatsMode::PpvNpv   => s.npv,
                }).collect();
                vals2.sort_by(|a, b| a.partial_cmp(b).unwrap());
                if !vals2.is_empty() {
                    let q1 = quantile(&vals2, 0.25);
                    let med = quantile(&vals2, 0.5);
                    let q3 = quantile(&vals2, 0.75);
                    let iqr = q3 - q1;
                    let lw = *vals2.iter()
                        .filter(|v| **v >= q1 - 1.5 * iqr)
                        .min_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap();
                    let uw = *vals2.iter()
                        .filter(|v| **v <= q3 + 1.5 * iqr)
                        .max_by(|a, b| a.partial_cmp(b).unwrap())
                        .unwrap();
                    let center_y2 = idx as f64 + 0.7;
                    let y0 = center_y2 - 0.1;
                    let y1 = center_y2 + 0.1;
                    chart.draw_series(std::iter::once(
                        Rectangle::new([(q1, y0), (q3, y1)], box_style.clone())
                    ))
                    .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
                    chart.draw_series(std::iter::once(
                        PathElement::new(vec![(med, y0), (med, y1)], whisker_style.clone())
                    ))
                    .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
                    chart.draw_series(std::iter::once(
                        PathElement::new(vec![(lw, center_y2), (q1, center_y2)], whisker_style.clone())
                    ))
                    .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
                    chart.draw_series(std::iter::once(
                        PathElement::new(vec![(q3, center_y2), (uw, center_y2)], whisker_style.clone())
                    ))
                    .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
                }
            }
        }
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


    let color_3 = match col3 {
        Some(color) => color.to_owned(),
        None => RGBColor(70, 70, 70)
    };

    let default_ref1 = data.values()
        .filter_map(|d| d.stats.iter()
            .find(|s| s.name == "consensus")
            .map(|s| s.sensitivity))
        .fold(0.0_f64, f64::max);

    let default_ref2 = data.values()
        .filter_map(|d| d.stats.iter()
            .find(|s| s.name == "consensus")
            .map(|s| s.specificity))
        .fold(0.0_f64, f64::max);

    let r1 = ref1.unwrap_or(default_ref1);
    let r2 = ref2.unwrap_or(default_ref2);

    chart.draw_series(DashedLineSeries::new(
        vec![(r1, 0.0), (r1, n_cats as f64)],
        5, 10,
        ShapeStyle {
            color: color_1.into(),
            filled: false,
            stroke_width: 2,
        },
    )).map_err(|e| CiqaError::StripPlotError(e.to_string()))?;


    // 5) draw dashed ref2 (max spec) line
    chart.draw_series(DashedLineSeries::new(
        vec![(r2, 0.0), (r2, n_cats as f64)],
        5, 10,
        ShapeStyle {
            color: color_2.into(),
            filled: false,
            stroke_width: 2,
        },
    )).map_err(|e| CiqaError::StripPlotError(e.to_string()))?;


    if ci {
        // === Draw the 95% CI Boxes behind the points ===

        let ci_box_style = ShapeStyle {
            color: RGBAColor(211, 211, 211, 0.3),
            stroke_width: 0,
            filled: true
        };

        for (idx, cat_name) in categories.iter().enumerate() {
            if let Some(diag_data) = data.get(cat_name) {
                let summary = &diag_data.summary;
                // Select the correct pair of CI intervals based on the mode.
                let (ci_m1, ci_m2) = match mode {
                    StatsMode::SensSpec => (summary.ci_sensitivity.clone(), summary.ci_specificity.clone()),
                    StatsMode::PpvNpv   => (summary.ci_ppv.clone(), summary.ci_npv.clone()),
                };

                // We span a vertical band from y = idx + 0.2 to idx + 0.4.
                chart.draw_series(std::iter::once(
                    Rectangle::new(
                        [(ci_m1[0], idx as f64 + 0.3), (ci_m1[1], idx as f64 + 0.5)],
                        ci_box_style.clone(),
                    ),
                ))
                .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;

                // We span a vertical band from y = idx + 0.6 to idx + 0.8.
                chart.draw_series(std::iter::once(
                    Rectangle::new(
                        [(ci_m2[0], idx as f64 + 0.5), (ci_m2[1], idx as f64 + 0.7)],
                        ci_box_style.clone(),
                    ),
                ))
                .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
            }
        }
    }

    if barplot {
        
        // First bar (e.g., Sensitivity/PPV): medium gray
        let bar_fill_1 = ShapeStyle {
            color: RGBColor(color_1.0, color_1.1, color_1.2).to_rgba().mix(0.6), // RGBAColor(100, 100, 100, 0.3),
            stroke_width: 0,
            filled: true,
        };

        // Second bar (e.g., Specificity/NPV): lighter gray
        let bar_fill_2 = ShapeStyle {
            color: RGBColor(color_2.0, color_2.1, color_2.2).to_rgba().mix(0.6), // RGBAColor(180, 180, 180, 0.3)
            stroke_width: 0,
            filled: true,
        };

        // Border color matches stripplot points
        let bar_border_1 = ShapeStyle {
            color: BLACK.to_rgba(), // RGBColor(color_1.0, color_1.1, color_1.2).to_rgba().mix(0.9),
            stroke_width: 1,
            filled: false,
        };
    
        let bar_border_2 = ShapeStyle {
            color: BLACK.to_rgba(),  // RGBColor(color_2.0, color_2.1, color_2.2).to_rgba().mix(0.9),
            stroke_width: 1,
            filled: false,
        };
    
        for (idx, cat_name) in categories.iter().enumerate() {
            if let Some(diag_data) = data.get(cat_name) {
                let stats = &diag_data.stats;
                let (mut vals1, mut vals2) = (Vec::new(), Vec::new());
                for s in stats {
                    match mode {
                        StatsMode::SensSpec => {
                            vals1.push(s.sensitivity);
                            vals2.push(s.specificity);
                        }
                        StatsMode::PpvNpv => {
                            vals1.push(s.ppv);
                            vals2.push(s.npv);
                        }
                    }
                }
    
                let avg1 = if !vals1.is_empty() {
                    vals1.iter().copied().sum::<f64>() / vals1.len() as f64
                } else { 0.0 };
    
                let avg2 = if !vals2.is_empty() {
                    vals2.iter().copied().sum::<f64>() / vals2.len() as f64
                } else { 0.0 };
    
                let bar_height = 0.06;
                let (bar1_y, bar2_y) = (idx as f64 + 0.4, idx as f64 + 0.6);
                let y0_1 = bar1_y - bar_height;
                let y1_1 = bar1_y + bar_height;
                let y0_2 = bar2_y - bar_height;
                let y1_2 = bar2_y + bar_height;
    
                chart.draw_series(vec![
                    Rectangle::new([(0.0, y0_1), (avg1, y1_1)], bar_fill_1.clone()),
                    Rectangle::new([(0.0, y0_1), (avg1, y1_1)], bar_border_1.clone()),
                    Rectangle::new([(0.0, y0_2), (avg2, y1_2)], bar_fill_2.clone()),
                    Rectangle::new([(0.0, y0_2), (avg2, y1_2)], bar_border_2.clone()),
                ])
                .map_err(|e| CiqaError::StripPlotError(e.to_string()))?;
            }
        }
    }

    // We'll use a small random jitter so multiple points don't overlap exactly
    let mut rng = rand::rng();

    let mut measure1_shapes = Vec::new();
    let mut measure2_shapes = Vec::new();

    let mut consensus_shapes = Vec::new();

    for (idx, cat_name) in categories.iter().enumerate() {

        let stats_vec = match data.get(cat_name) {
            Some(v) => v.stats.clone(),
            None => continue,
        };

        for stat in stats_vec {


            let (x1, x2) = match mode {
                StatsMode::SensSpec => (stat.sensitivity, stat.specificity),
                StatsMode::PpvNpv   => (stat.ppv, stat.npv),
            };

            // y offsets for measure1 vs measure2
            let base_1 = idx as f64 + 0.4;
            let base_2 = idx as f64 + 0.6;

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
        .position(legend_position.unwrap_or(SeriesLabelPosition::MiddleMiddle))
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


pub enum CellShape {
    Circle,
    Square { border_width: u32 },
}

pub fn plot_diagnostic_matrix(
    data: &[Vec<DiagnosticReview>],
    reference: Option<&[DiagnosticReview]>,
    consensus: Option<&[DiagnosticReview]>,
    output: &Path,
    palette: &Palette,
    shape: CellShape,
    column_header: PanelColumnHeader,
    title: Option<&str>,
    header_text: Option<&str>,
    consensus_stats: Option<&DiagnosticStats>
) -> Result<(), CiqaError> {
    // 1. Determine sample IDs and sanity‑check all columns
    let sample_labels: Vec<_> = if let Some(col) = reference {
        col.iter().map(|r| r.sample_id.clone()).collect()
    } else if !data.is_empty() {
        data[0].iter().map(|r| r.sample_id.clone()).collect()
    } else if let Some(col) = consensus {
        col.iter().map(|r| r.sample_id.clone()).collect()
    } else {
        return Err(CiqaError::NoColumnsFound);
    };
    let nrows = sample_labels.len();
    for col in data {
        assert_eq!(col.len(), nrows, "Data column length mismatch");
    }
    if let Some(col) = reference {
        assert_eq!(col.len(), nrows, "Reference length mismatch");
    }
    if let Some(col) = consensus {
        assert_eq!(col.len(), nrows, "Consensus length mismatch");
    }

    // 2. Flatten into a single Vec<&[DiagnosticReview]> in draw order
    let mut columns: Vec<&[DiagnosticReview]> = Vec::new();
    for col in data { columns.push(col.as_slice()); }
    if let Some(col) = consensus { columns.push(col); }
    if let Some(col) = reference { columns.push(col); }

    // meta column index (always last)
    let meta_idx = columns.len();
    let total_cols = columns.len() + 1;

    // 3. Panel layout parameters (same as QC function)
    let chunk_size   = 12;
    let max_panels_x = 4;
    let num_panels   = (nrows + chunk_size - 1) / chunk_size;
    let panels_x     = std::cmp::min(max_panels_x, num_panels);
    let panels_y     = (num_panels + max_panels_x - 1) / max_panels_x;

    // 4. Margins, paddings & cell sizing
    let outer_margin        = 20.0;
    let panel_padding_left  = 210.0;
    let panel_padding_right = 40.0;
    let panel_padding_bottom= 10.0;
    let col_padding         = 4.0;
    let row_padding         = 4.0;
    let cell_px             = 24;
    let font_size           = (cell_px as f64).clamp(4.0, 14.0).round() as u32;
    let legend_height_px    = 120.0;

    // 5. Compute top‑paddings per panel (for column headers)
    let base_padding_y   = 10;
    let header_padding_y = 100;
    let (panel_paddings_top, total_padding_top) =
        column_header.panel_paddings_top(base_padding_y, header_padding_y, panels_x, panels_y);
    let panel_column_headers =
        column_header.panel_column_header(panels_x, panels_y);

    // 6. Figure out total SVG size

    let stride = cell_px as f64 + col_padding;
    let consensus_gap_px = stride; // one-cell gap before Consensus
    
    let consensus_idx: Option<usize> = if consensus.is_some() { Some(data.len()) } else { None };
    let extra_gap_per_panel = if consensus_idx.is_some() { consensus_gap_px } else { 0.0 };
    
    let height_px = (
        (panels_y as f64 * chunk_size as f64) * (cell_px as f64 + row_padding)
        + (panels_y as f64 * panel_padding_bottom)
        + total_padding_top as f64
        + legend_height_px                         // <— reserve legend space
        + (2.0 * outer_margin)
    ).ceil() as u32;
    
    let width_px = (
        (panels_x as f64 * total_cols as f64) * stride
        + (panels_x as f64 * extra_gap_per_panel)
        + (panels_x as f64 * (panel_padding_left + panel_padding_right))
        + (2.0 * outer_margin)
    ).ceil() as u32;

    // 7. Start drawing
    let root = SVGBackend::new(output, (width_px, height_px)).into_drawing_area();
    root.fill(&WHITE)?;

    // 8. Optional title
    if let Some(t) = title {
        root.draw_text(
            t,
            &("monospace", 16).into_font().into_text_style(&root)
                .pos(Pos::new(HPos::Center, VPos::Top)),
            ((width_px / 2) as i32, 20),
        )?;
    }

    // 9. Split into panels
    let plot_area = root.margin(outer_margin, outer_margin, outer_margin, outer_margin);
    let panels = plot_area.split_evenly((panels_y, panels_x));

    // 10. Style lookup for each review
    let get_style = |rev: &DiagnosticReview| {
        palette.colors[rev.outcome.index()].filled()
    };

    // helper: x position with gap and meta
    let x_for = |col_idx: usize| -> f64 {
        let base = (col_idx as f64) * stride;
        match consensus_idx {
            Some(ci) if col_idx >= ci => base + consensus_gap_px,
            _ => base,
        }
    };

    // 11. Draw each panel
    for (panel, (chunk_idx, chunk)) in
        panels.iter().zip(sample_labels.chunks(chunk_size).enumerate())
    {
        let top_pad = panel_paddings_top[chunk_idx] as f64;
        let panel_area = panel.margin(
            top_pad,
            panel_padding_bottom,
            panel_padding_left,
            panel_padding_right,
        );

        let nrows_chunk = chunk.len();

        // 11a. Draw every column (data → consensus → reference)
        for (col_idx, col) in columns.iter().enumerate() {
            let x0 = x_for(col_idx);
            for row in 0..nrows_chunk {
                let row_idx = chunk_idx * chunk_size + row;
                if row_idx >= nrows { continue; }
                let rev = &col[row_idx];
                let style = get_style(rev);
                let y0 = row as f64 * (cell_px as f64 + row_padding);
                match shape {
                    CellShape::Circle => {
                        let cx = (x0 + cell_px as f64 / 2.0) as i32;
                        let cy = (y0 + cell_px as f64 / 2.0) as i32;
                        panel_area.draw(&Circle::new((cx, cy), (cell_px as f64 / 2.0) as i32, style.clone()))?;
                    }
                    CellShape::Square { border_width } => {
                        let rect = [(x0 as i32, y0 as i32), ((x0 + cell_px as f64) as i32, (y0 + cell_px as f64) as i32)];
                        panel_area.draw(&Rectangle::new(rect, style.clone()))?;
                        panel_area.draw(&Rectangle::new(rect, WHITE.stroke_width(border_width)))?;
                    }
                }
            }
        }

        // 11b. Sample‐ID labels
        for (row, label) in chunk.iter().enumerate() {
            let y0 = row as f64 * (cell_px as f64 + row_padding) + (cell_px as f64 / 2.0) + row_padding;
            panel_area.draw_text(
                label,
                &("monospace", font_size)
                    .into_font()
                    .into_text_style(&panel_area)
                    .pos(Pos::new(HPos::Right, VPos::Center)),
                (-5, y0 as i32),
            )?;
        }

        // 11c. Rotated column headers on the first panel row
        if panel_column_headers[chunk_idx] {
            for col_idx in 0..columns.len() {           // only real columns get rotated headers
                let x0 = x_for(col_idx);
                let header = if col_idx < data.len() {
                    format!("{} {}", header_text.unwrap_or("Replicate"), col_idx + 1)
                } else if col_idx < data.len() + consensus.as_ref().map(|_| 1).unwrap_or(0) {
                    "Consensus".into()
                } else {
                    "Reference".into()
                };
                panel_area.draw_text(
                    &header,
                    &("monospace", font_size)
                        .into_font()
                        .into_text_style(&panel_area)
                        .pos(Pos::new(HPos::Left, VPos::Center))
                        .transform(FontTransform::Rotate270),
                    ((x0 + (cell_px as f64 / 2.0)) as i32, -5),
                )?;
            }

        }
    }
        // styles once
    let normal = ("monospace", font_size).into_font();
    let bold   = ("monospace", font_size).into_font().style(FontStyle::Bold);
    let gap_small  = 10i32;
    // helper: pixel width of text with a style
    let text_w = |s: &str, bolding: bool| -> i32 {
        let ts = if bolding { bold.clone() } else { normal.clone() }
            .into_text_style(&root);
        root.estimate_text_size(s, &ts).unwrap_or((((font_size as f64) * 0.60).round() as u32, 0.0 as u32)).0 as i32
    };

    // ---- measure total width precisely (center the whole line) ----
    let sw = 14.0;

    let legend = [
            ("TP", DiagnosticOutcome::TruePositive),
            ("TN", DiagnosticOutcome::TrueNegative),
            ("FP", DiagnosticOutcome::FalsePositive),
            ("FN", DiagnosticOutcome::FalseNegative),
            ("Control", DiagnosticOutcome::Control),
            ("Failed", DiagnosticOutcome::Indeterminate),
            ("N/A", DiagnosticOutcome::NotConsidered)
    ];

    let mut total_w: i32 = 0;
    for (i, (label, _)) in legend.iter().enumerate() {
        total_w += (sw as i32) + gap_small + text_w(label, false);
        if i != legend.len() - 1 { total_w += gap_small; }
    }

    let mut stats_fmt: Option<(String,String,String,String,String,String)> = None;
    if let Some(cs) = consensus_stats {
        let sens = format!("{:.1}%", cs.sensitivity * 100.0);
        let spec = format!("{:.1}%", cs.specificity * 100.0);
        let ppv  = format!("{:.1}%", cs.ppv * 100.0);
        let npv  = format!("{:.1}%", cs.npv * 100.0);
        let rc: String = format!("{:.1}%", average_replicate_certainty(data, reference));

        let nstr = format!("n = {}", cs.total);

        let pairs = [
            format!("Sensitivity: {sens}"),
            format!("Specificity: {spec}"),
            format!("PPV: {ppv}"),
            format!("NPV: {npv}"),
            format!("Replicate Certainty: {rc}")
        ];

        for (i, val) in pairs.iter().enumerate() {
            total_w += text_w(val, false) + gap_small;
        }

        total_w += text_w(&nstr, false) + gap_small;
        
        stats_fmt = Some((sens, spec, ppv, npv, nstr, rc));
    }

    let baseline_y = (height_px as f64 - outer_margin - 8.0).round() as i32;
    let mut x = ((width_px as i32) - total_w) / 2;

    // ---- draw: icons + labels ----
    for (label, out) in legend {
        let style = palette.colors[out.index()].filled();
        let icon_top = baseline_y as f64 - sw;

        match shape {
            CellShape::Circle => {
                let cx = (x as f64 + sw/2.0) as i32;
                let cy = (icon_top + sw/2.0) as i32;
                root.draw(&Circle::new((cx, cy), (sw/2.0) as i32, style.clone()))?;
            }
            CellShape::Square { .. } => {
                let r = [(x, icon_top as i32), ((x as f64 + sw) as i32, baseline_y)];
                root.draw(&Rectangle::new(r, style.clone()))?;
            }
        }
        x += (sw as i32) + gap_small;

        root.draw_text(
            label,
            &normal.clone().into_text_style(&root).pos(Pos::new(HPos::Left, VPos::Center)),
            (x, baseline_y),
        )?;
        x += text_w(label, false) + gap_small;
    }

    // ---- draw: stats ----
    if let Some((sens, spec, ppv, npv, nstr, rc)) = stats_fmt {
        let mut draw_pair = |val: &str| -> Result<(), CiqaError> {
            root.draw_text(
                val,
                &normal.clone().into_text_style(&root).pos(Pos::new(HPos::Left, VPos::Center)),
                (x, baseline_y),
            )?;
            x += text_w(val, false) + gap_small;
            Ok(())
        };

        draw_pair(&format!("Sensitivity: {sens}"))?;
        draw_pair(&format!("Specificity: {spec}"))?;
        draw_pair(&format!("PPV: {ppv}"))?;
        draw_pair(&format!("NPV: {npv}"))?;
        draw_pair(&format!("Replicate Certainty: {rc}"))?;
        draw_pair(&nstr)?;

    }

    Ok(())
}


#[derive(Serialize, Deserialize, Debug, PartialEq, Clone, clap::ValueEnum)]
pub enum PanelColumnHeader {
    Panel,
    FirstRow,
}
impl PanelColumnHeader {
    pub fn panel_paddings_top(
        &self,
        base_padding: usize,
        header_padding: usize,
        num_panels_x: usize,
        num_panels_y: usize,
    ) -> (Vec<usize>, usize) {
        let num_panels = num_panels_x * num_panels_y;

        let paddings: Vec<usize> = match self {
            PanelColumnHeader::Panel => {
                // every panel gets the header+base
                vec![base_padding + header_padding; num_panels]
            }
            PanelColumnHeader::FirstRow => {
                // only the *first row* of panels
                let mut v = Vec::with_capacity(num_panels);
                // panels_x panels in row 0 get header
                v.extend(std::iter::repeat(base_padding + header_padding).take(num_panels_x));
                // the rest get just base
                v.extend(std::iter::repeat(base_padding).take(num_panels - num_panels_x));
                v
            }
        };

        let total_padding_y = match self {
            PanelColumnHeader::Panel => (base_padding+header_padding)*num_panels_y,
            PanelColumnHeader::FirstRow => (base_padding+header_padding)+(base_padding*(num_panels_y-1))
        };

        (paddings, total_padding_y)
    }
    pub fn panel_column_header(&self, num_panels_x: usize, num_panels_y: usize) -> Vec<bool>{
        let num_panels = num_panels_x*num_panels_y;
        match self {
            PanelColumnHeader::Panel => vec![true; num_panels],
            PanelColumnHeader::FirstRow => {
                [vec![true; num_panels_x], vec![false; num_panels-num_panels_x]].concat()
            }
        }
    }
}


pub fn plot_qc_summary_matrix(
    summaries: &[QualityControlSummary],
    overall: Option<&[QcStatus]>,
    output: &Path,
    shape: CellShape,
    column_header: PanelColumnHeader,
    title: Option<&str>,
) -> Result<(), CiqaError> {

    let mut sorted_summaries = summaries.to_vec();
    sorted_summaries.sort_by_key(|s| s.id.clone());

    let sample_labels: Vec<_> = sorted_summaries.iter().map(|s| s.id.clone()).collect();

    let mut categories = vec![
        "Input Reads",
        "Output Reads",
        "ERCC | EDCC",
        "Read Quality",
    ];

    if sorted_summaries.iter().any(|s| s.phage_coverage.is_some()) {
        categories.push("Phage Coverage");
    }

    let has_overall = overall.is_some();
    let total_cols = categories.len() + if has_overall { 1 } else { 0 };


    let nrows = sample_labels.len();
    let chunk_size = 12;
    let max_panels_x = 3;

    let num_panels = (nrows + chunk_size - 1) / chunk_size;
    let panels_x = std::cmp::min(max_panels_x as usize, num_panels);
    let panels_y = (num_panels + max_panels_x as usize - 1) / max_panels_x as usize;

    
    let outer_margin = 20.0;         
    let panel_padding_left = 210.0; 
    let panel_padding_right = 40.0;
    let panel_padding_bottom = 10.0;   
    
    let col_padding = 4.0;    
    let row_padding = 4.0;    
    let cell_px = 24;                                                                   // target pixel height/width per row
    let font_size = (cell_px as f64).clamp(4.0, 14.0).round() as u32;

    
    let total_padding_top = match column_header {
        PanelColumnHeader::FirstRow => 0,  // added to margin
        PanelColumnHeader::Panel => 100*panels_y
    };

    let panel_column_headers =
        column_header.panel_column_header(panels_x, panels_y);

    let height_px = (
        (panels_y as f64 * chunk_size as f64) * (cell_px as f64 + row_padding)
        + (panels_y as f64 * panel_padding_bottom)
        + total_padding_top as f64
        + (2.0 * outer_margin)
    ).ceil() as u32;

    let width_px = (
        (panels_x as f64 * total_cols as f64) * (cell_px as f64 + col_padding) 
        + (panels_x as f64 * (panel_padding_left+panel_padding_right)) 
        + (2.0 * outer_margin)
    ).ceil() as u32;

    log::info!("Total: {nrows} Panels X: {panels_x} Panels Y: {panels_y} Outer margin: {outer_margin} Inter-column padding: {col_padding} Cell size in px: {cell_px} Height in px: {height_px}  Width in px: {width_px}");

    let root = SVGBackend::new(
        output, 
        (width_px, height_px as u32)
    ).into_drawing_area();

    root.fill(&WHITE)?;

    if let Some(t) = title {
        root.draw_text(
            t,
            &("monospace", 16).into_font().into_text_style(&root)
                .pos(Pos::new(HPos::Center, VPos::Top)),
            ((width_px / 2) as i32, 20),
        )?;
    }

    let outer_margin_top = match column_header {
        PanelColumnHeader::FirstRow => outer_margin+100.0,  // added to margin
        PanelColumnHeader::Panel => outer_margin
    };

    let plot_area = root.margin(outer_margin_top, outer_margin, outer_margin, outer_margin);
    let panels = plot_area.split_evenly((panels_y, panels_x));

    let get_style = |status: &QcStatus| match status {
        QcStatus::Fail => RGBColor(204,  32,  32).filled(),
        QcStatus::Pass => RGBColor(  0, 102, 204).filled(),
        QcStatus::Ok   => RGBColor(  0, 153,   0).filled(),
    };

    for (panel, (chunk_idx, chunk)) in panels.iter().zip(sample_labels.chunks(chunk_size).enumerate()) {
        
        let panel_padding_top = match column_header {
            PanelColumnHeader::FirstRow => 0,  // added to margin
            PanelColumnHeader::Panel => 100
        };
        
        log::info!("Panel {chunk_idx} padding top: {panel_padding_top} with total padding (y) {total_padding_top}");

        let panel_area = panel.margin(
            panel_padding_top as f64, 
            panel_padding_bottom, 
            panel_padding_left, 
            panel_padding_right
        );

        let nrows_chunk = chunk.len();
        
        let stride   = cell_px as f64 + col_padding;

        for (col_idx, category) in categories.iter().enumerate() {

            let x0 = (col_idx as f64) * stride;
            for row in 0..nrows_chunk {
                let row_idx = chunk_idx * chunk_size + row;
                if row_idx >= sorted_summaries.len() {
                    continue;
                }
                let summary = &sorted_summaries[row_idx];

                let status = match *category {
                    "Input Reads"     => &summary.input_reads,
                    "Output Reads"    => &summary.output_reads,
                    "ERCC | EDCC"     => &summary.ercc_constructs,
                    "Read Quality"    => &summary.fastp_status,
                    "Phage Coverage"  => summary.phage_coverage.as_ref().unwrap(),
                    _ => continue,
                };
                let style = get_style(status);

                let y0 = row as f64 * (cell_px as f64 + row_padding);
                match shape {
                    CellShape::Circle => {
                        let cx = (x0 + cell_px as f64 / 2.0) as i32;
                        let cy = (y0 + cell_px as f64 / 2.0) as i32;
                        panel_area.draw(&Circle::new((cx, cy), (cell_px as f64 / 2.0) as i32, style.clone()))?;
                    }
                    CellShape::Square { border_width } => {
                        let rect = [
                            ((x0) as i32, y0 as i32),
                            ((x0 + cell_px as f64) as i32, (y0 + cell_px as f64) as i32),
                        ];
                        panel_area.draw(&Rectangle::new(rect, style.clone()))?;
                        panel_area.draw(&Rectangle::new(rect, WHITE.stroke_width(border_width)))?;
                    }
                }
            }
        }
        // Sample‐ID labels on the left edge of the panel
        for (row, label) in chunk.iter().enumerate() {
            let y0 = row as f64 * (cell_px as f64 + row_padding) + (cell_px as f64 / 2.0) + row_padding;
            panel_area.draw_text(
                label,
                &("monospace", font_size)
                    .into_font()
                    .into_text_style(&panel_area)
                    .pos(Pos::new(HPos::Right, VPos::Center)),
                (-5, y0 as i32),
            )?;
        }
        if panel_column_headers[chunk_idx] {
            for (col_idx, category) in categories.iter().enumerate() {
                let x0 = (col_idx as f64) * stride;
                panel_area.draw_text(
                    category,
                    &("monospace", font_size)
                        .into_font()
                        .into_text_style(&panel_area)
                        .pos(Pos::new(HPos::Left, VPos::Bottom))
                        .transform(FontTransform::Rotate270),
                    ((x0 + (cell_px as f64 / 2.0)) as i32, -5), // y = 0 is top of panel_area due to margin
                )?;
            }
        }
        
        // Overall column
        if has_overall {

            let overall_idx = categories.len();
            let x0: f64 = (overall_idx as f64) * stride;

            for row in 0..nrows_chunk {
                let row_idx = chunk_idx * chunk_size + row;
                if row_idx >= sorted_summaries.len() {
                    continue;
                }
                let status = &overall.unwrap()[row_idx];
                let style = get_style(status);

                let y0 = row as f64 * (cell_px as f64 + row_padding);
                match shape {
                    CellShape::Circle => {
                        let cx = (x0 + cell_px as f64 / 2.0) as i32;
                        let cy = (y0 + cell_px as f64 / 2.0) as i32;
                        panel_area.draw(&Circle::new((cx, cy), (cell_px as f64 / 2.0) as i32, style.clone()))?;
                    }
                    CellShape::Square { border_width } => {
                        let rect = [
                            ((x0) as i32, y0 as i32),
                            ((x0 + cell_px as f64) as i32, (y0 + cell_px as f64) as i32),
                        ];
                        panel_area.draw(&Rectangle::new(rect, style.clone()))?;
                        panel_area.draw(&Rectangle::new(rect, WHITE.stroke_width(border_width)))?;
                    }
                }
            }
        }

        
    }

    Ok(())
}