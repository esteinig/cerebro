use std::path::PathBuf;

#[cfg(feature = "local")]
use meta_gpt::{gpt::TreeConfig, model::GeneratorModel};
use meta_gpt::gpt::{AssayContext, AgentPrimer, TaskConfig};

use clap::{Args, Parser, Subcommand};
use cerebro_model::api::cerebro::schema::SampleType;
use crate::plate::{PanelColumnHeader, MissingOrthogonal, StatsMode};

/// Cerebro: production file system watcher 
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-ciqa", version)]
pub struct App {
    /// API URL
    #[clap(
        long, 
        short = 'u', 
        default_value = "http://api.cerebro.localhost", 
        env = "CEREBRO_API_URL"
    )]
    pub url: String,
    /// API token - usually provided with CEREBRO_API_TOKEN
    #[clap(
        long, 
        short = 'e', 
        env = "CEREBRO_API_TOKEN",
        hide_env_values = true
    )]
    pub token: Option<String>,
    /// API token file - can be set from environment variable
    #[clap(
        long, 
        short = 'f', 
        env = "CEREBRO_API_TOKEN_FILE"
    )]
    pub token_file: Option<PathBuf>,
    /// User team name or identifier for requests that require team specification 
    #[clap(
        long, 
        short = 't', 
        env = "CEREBRO_USER_TEAM",
        hide_env_values = true
    )]
    pub team: Option<String>,
    /// Team database name or identifier for requests that require database access 
    #[clap(
        long, 
        short = 'd', 
        env = "CEREBRO_USER_DB",
        hide_env_values = true
    )]
    pub db: Option<String>,
    /// Team database project name or identifier for requests that require project access 
    #[clap(
        long, 
        short = 'p', 
        env = "CEREBRO_USER_PROJECT",
        hide_env_values = true
    )]
    pub project: Option<String>,
    /// SeaweedFS master node address
    #[clap(
        long, 
        short = 'a',
        default_value = "http://fs.cerebro.localhost", 
        env = "CEREBRO_FS_URL"
    )]
    pub fs_url: String,
    /// SeaweedFS master node port
    #[clap(
        long, 
        short = 'm',
        env = "CEREBRO_FS_PORT",
        default_value = "9333", 
    )]
    pub fs_port: String,
    /// SSL certificate verification is ignored [DANGER]
    #[clap(
        long, 
        env = "CEREBRO_DANGER_ACCEPT_INVALID_TLS_CERTIFICATE"
    )]
    pub danger_invalid_certificate: bool,
    
    #[clap(subcommand)]
    pub command: Commands,
}


#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Open TUI for inspecting output from reference plate evaluations
    Tui(TuiArgs),
    /// Prefetch tiered data with filter settings for a sample and save to file
    Prefetch(PrefetchArgs),
    /// Write the plate to a tab-delimited table 
    WriteTable(WriteTableArgs),
    /// Plot the tiered prefetch taxa composition
    PlotPrefetch(PlotPrefetchArgs),
    /// Upload a CI/QA dataset to Cerebro Fs and index/stage for pipeline
    UploadCiqaDataset(UploadCiqaDatasetArgs),
    /// Plot the reference plate layout
    PlotPlate(PlotPlateArgs),
    /// Plot the reference plate quality control summary
    PlotQc(PlotQcArgs),
    /// Plot the reference plate review results
    PlotReview(PlotReviewArgs),
    /// Review diagnostic outcome of a review table against the reference data
    Review(ReviewArgs),
    /// Compare two diagnostic consensus reviews using McNemar's test statistic
    Mcnemar(McnemarArgs),
    #[cfg(feature = "local")]
    /// Diagnose samples on the reference plate using the generative practitioner
    DiagnoseLocal(DiagnoseLocalArgs),
    /// Debug the pathogen calls made in a set of diagnostic practitioner outputs
    DebugPathogen(DebugPathogenArgs),
}

#[derive(Debug, Args)]
pub struct WriteTableArgs {
    /// Reference plate file (.json)
    #[clap(long, short = 'p')]
    pub plate: PathBuf,
    /// Output plot file (plate.tsv)
    #[clap(long, short = 'o', default_value="plate.tsv")]
    pub output: PathBuf,
    /// Include only positive taxa if test included detection at species rank
    #[clap(long, short = 's')]
    pub species_rank: bool,
}


#[derive(Debug, Args)]
pub struct McnemarArgs {
    /// Path to output from review of configuration A (DiagnosticData, JSON)
    #[clap(long, short = 'a')]
    pub review_a: String,
    /// Path to output from review of configuration B (DiagnosticData, JSON)
    #[clap(long, short = 'b')]
    pub review_b: String,
}

#[derive(Debug, Args)]
pub struct PlotPlateArgs {
    /// File paths for 
    #[clap(long, short = 's')]
    pub sample: String,
}


#[derive(Debug, Args)]
pub struct PlotQcArgs {
    /// Cerebro model files from DB (.json)
    #[clap(long, short = 'i', num_args=0..)]
    pub cerebro: Vec<PathBuf>,
    /// Output plot file (.svg)
    #[clap(long, short = 'o', default_value="qc_summary.svg")]
    pub output: PathBuf,
    /// Output positive controls summary file (.tsv)
    #[clap(long, short = 'p')]
    pub positive_controls: Option<PathBuf>,
    /// Output directory for summary files (.json)
    #[clap(long, short = 'd',)]
    pub outdir: Option<PathBuf>,
    /// Threads to use for reading models
    #[clap(long, short = 't', default_value="4")]
    pub threads: u64,
    /// Column header format
    #[clap(long, default_value="panel")]
    pub column_header: PanelColumnHeader,
    /// Set plot title
    #[clap(long)]
    pub title: Option<String>,
}


#[derive(Debug, Args)]
pub struct PlotPrefetchArgs {
    /// Prefetch file (.json)
    #[clap(long, short = 'i')]
    pub prefetch: PathBuf,
    /// Output plot file (.svg)
    #[clap(long, short = 'o', default_value="qc_summary.svg")]
    pub output: PathBuf,
    /// Plot width (px)
    #[clap(long, default_value="800")]
    pub width: u32,
    /// Plot height (px)
    #[clap(long, default_value="600")]
    pub height: u32,
    /// Set plot title
    #[clap(long)]
    pub title: Option<String>,
}

#[derive(Debug, Args)]
pub struct UploadCiqaDatasetArgs {
    
    /// Read files dataset with NTC/ENV control for showread CI/QA
    #[clap(long, short = 'f', num_args=1..)]
    pub fastq_pe: Vec<PathBuf>,
    /// Sequence run identifier
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Biological sample identifier
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,
    /// Pipeline run identifier
    #[clap(long, short = 'p')]
    pub pipeline_id: Option<String>,
    /// File description
    #[clap(long, short = 'd')]
    pub description: Option<String>,
    
}

#[derive(Debug, Args)]
pub struct DebugPathogenArgs {
    /// Diagnostic agent output files (.json)
    #[clap(long, short = 'g', num_args=1..)]
    pub gpt: Vec<PathBuf>,
    /// Summary of pathogen calls output (.tsv)
    #[clap(long, short = 'o')]
    pub output: PathBuf
}

#[derive(Debug, Args)]
pub struct PlotReviewArgs {
    /// Review stats (.json) for a category - file stem is the category name
    #[clap(long, short = 's', num_args=1..)]
    pub stats: Vec<PathBuf>,
    /// Diagnostic statistics selection
    #[clap(long, short = 'm', default_value="sens-spec")]
    pub mode: StatsMode,
    /// Output plot file (.svg)
    #[clap(long, short = 'o', default_value="review_sens_spec.svg")]
    pub output: PathBuf,
    /// Plot 95% CI interval as grey box behind points (assuming normal distribution)
    #[clap(long, short = 'c')]
    pub ci: bool,
    /// Plot boxplot overlay
    #[clap(long, short = 'b')]
    pub boxplot: bool,
    /// Plot barplot overlay
    #[clap(long, short = 'a')]
    pub barplot: bool,
    /// Custom y-labels 
    #[clap(long, short = 'y', num_args=0..)]
    pub y_labels: Option<Vec<String>>,
    /// Reference line 1 (sensitivity/ppv)
    #[clap(long)]
    pub ref1: Option<f64>,
    /// Reference line 1 (specificty/npv)
    #[clap(long)]
    pub ref2: Option<f64>,
    /// Plot width (px)
    #[clap(long, default_value="800")]
    pub width: u32,
    /// Plot height (px)
    #[clap(long, default_value="600")]
    pub height: u32,
}


#[cfg(feature = "local")]
#[derive(Debug, Args)]
pub struct DiagnoseLocalArgs {
    /// Reference plate file (.json)
    #[clap(long, short = 'p')]
    pub plate: PathBuf,
    /// Data directory for prefetched data files (.prefetch.json)
    #[clap(long, short = 'd')]
    pub prefetch: PathBuf,
    /// Output directory for diagnostic results (.json)
    #[clap(long, short = 'o')]
    pub outdir: PathBuf,
    /// Local generator model for diagnostic queries
    #[clap(long, short = 'm', default_value="deepseekr1-qwen7b-q8-0")]
    pub model: GeneratorModel,
    /// The length of the sample to generate (in tokens) - needs to be long enough to capture model thought process
    #[arg(short = 'n', long, default_value_t = 20000)]
    pub sample_len: usize,
    /// The temperature used to generate samples, use 0 for greedy sampling
    #[arg(long, short='t', default_value_t = 0.8)]
    pub temperature: f64,
    /// Nucleus sampling probability cutoff.
    #[arg(long, short='s')]
    pub top_p: Option<f64>,
    /// Only sample among the top K samples.
    #[arg(long, short='k')]
    pub top_k: Option<usize>,
    /// Minimum probability threshold for sampling (min_p sampling).
    #[arg(long)]
    pub min_p: Option<f64>,
    /// Number of GPUs to batch sample evaluations on (each must support model)
    #[arg(long, short='g', default_value_t=1)]
    pub num_gpu: usize,
    /// Include clinical notes from plate reference into prompt context (if available)
    #[clap(long, short = 'c')]
    pub clinical_notes: bool,
    /// Include sample description from plate reference into prompt context
    #[clap(long, short = 's', default_value="true")]
    pub sample_context: Option<bool>,
    /// Include assay context into prompt context
    #[clap(long, short = 'a', default_value="cerebro-filter")]
    pub assay_context: Option<AssayContext>,
    /// System prompt primer for diagnostic agent
    #[clap(long, default_value="default")]
    pub agent_primer: Option<AgentPrimer>,
    /// Task configuration sets for the tiered decision tree
    #[clap(long, default_value="default")]
    pub task_config: TaskConfig,
    /// Task configuration sets for the tiered decision tree
    #[clap(long, default_value="tiered")]
    pub tree_config: TreeConfig,
    /// Post process taxa after filtering and retrieval by collapsing species variants and selecting best species per genus (Archaea|Bacteria|Eukaryota)
    #[clap(long, default_value="true")]
    pub post_filter: Option<bool>,
    /// Minimum species per genus required to enable selecting best species for the genus during post filter processing (Archaea|Bacteria|Eukaryota)
    #[clap(long, default_value="3")]
    pub min_species: usize,
    /// Override default species reduction filter to these domains (Archaea|Bacteria|Eukaryota)
    #[clap(long, num_args=1..)]
    pub species_domains: Option<Vec<String>>,
    /// Override default variant species collapse using pruned species name (GTDB, Archaea|Bacteria)
    #[clap(long, default_value="true")]
    pub collapse_variants: Option<bool>,
    /// Override default phage exclusion option for post processing taxa after filtering and retrieval
    #[clap(long, default_value="true")]
    pub exclude_phage: Option<bool>,
    /// Force overwrite output, otherwise skip if exists
    #[clap(long, short = 'f')]
    pub force: bool,
    /// Model directory for generator model(.gguf) and tokenizer file (.json)
    #[clap(long, default_value=".")]
    pub model_dir: PathBuf,
    /// Disable thinking in Qwen3
    #[clap(long)]
    pub disable_thinking: bool,
}


#[derive(Debug, Args)]
pub struct ReviewArgs {
    /// Reference plate schema (.json)
    #[clap(long, short = 'p')]
    pub plate: PathBuf,
    /// Plate review files (.tsv) or `meta-gpt` output directories with diagnostic result files if `--diagnostic_agent` flag enabled
    #[clap(long, short = 'r', num_args=1..)]
    pub review: Vec<PathBuf>,
    /// Output data (.json) from the reviews 
    /// 
    /// Stores an array of `DiagnosticStats` JSON schemas with attribute `name` assigned the review directory name
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Enable diagnostic meta-gpt agent output review
    #[clap(long, short = 'd')]
    pub diagnostic_agent: bool,
    /// Evaluate reference samples by sample type
    #[clap(long, short = 's')]
    pub sample_type: Option<SampleType>,
    /// Set reference results for these samples to 'None'
    #[clap(long, short = 'n', num_args=1..)]
    pub set_none: Option<Vec<String>>,
    /// Handle references with missing orthogonal data and a positive review call
    #[clap(long, short = 'm', default_value="indeterminate")]
    pub missing_orthogonal: MissingOrthogonal,
    /// Reference plate review (.json) for plot
    #[clap(long)]
    pub reference: Option<PathBuf>,
    /// Set plot width
    #[clap(long)]
    pub plot: Option<PathBuf>,
    /// Set plot width
    #[clap(long, default_value="950")]
    pub width: u32,
    /// Set plot height
    #[clap(long, default_value="600")]
    pub height: u32,
    /// Set plot title
    #[clap(long)]
    pub title: Option<String>,
}


#[derive(Debug, Args, Clone)]
pub struct PrefetchArgs {
    /// Reference plate file (.json)
    #[clap(long, short = 'p')]
    pub plate: PathBuf,
    /// Output directory for tiered filter data (.prefetch.json)
    #[clap(long, short = 'o')]
    pub outdir: PathBuf,
    /// Subset of sample identifiers from plate to prefetch only
    #[clap(long, short = 's', num_args(1..))]
    pub samples: Option<Vec<String>>,
    /// Prevalence contamination config (.json) otherwise default
    #[clap(long)]
    pub contamination: Option<PathBuf>,
    /// Tiered filter config (.json) otherwise default
    #[clap(long)]
    pub tiered_filter: Option<PathBuf>,
    /// Disable any filters or contamination controls, simply prefetch 
    /// the raw data into the primary filter category
    #[clap(long)]
    pub disable_filter: bool,
    /// Threads to use for fetching data
    #[clap(long, short = 't', default_value="4")]
    pub threads: u64,
    /// Force overwrite output, otherwise skip if exists
    #[clap(long, short = 'f')]
    pub force: bool,
    /// Write JSON summary of all prefetched samples
    #[clap(long)]
    pub summary: Option<PathBuf>,
    /// Write TSV of positives with no detected candidates
    #[clap(long)]
    pub summary_missed: Option<PathBuf>,
    /// Path to an existing unfiltered PathogenDetectionTable CSV to subset
    #[clap(long)]
    pub pathogen_detection_table: Option<PathBuf>,
    /// Output TSV for subset of PD table containing only missed detections
    #[clap(long)]
    pub pathogen_missed: Option<PathBuf>,
}


#[derive(Debug, Args, Clone)]
pub struct TuiArgs {
}

#[derive(Debug, Args)]
pub struct GlobalOptions {
    
}


pub fn get_styles() -> clap::builder::Styles {
	clap::builder::Styles::styled()
		.header(
			anstyle::Style::new()
				.bold()
				.underline()
				.fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
		)
		.literal(
			anstyle::Style::new()
				.bold()
				.fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
		)
}
