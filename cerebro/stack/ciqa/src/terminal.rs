use std::path::PathBuf;
use cerebro_gp::{gpt::{AssayContext, GptModel}, text::GeneratorModel};
use clap::{ArgGroup, Args, Parser, Subcommand};

use crate::plate::{DiagnosticOutcome, MissingOrthogonal, SampleType, StatsMode};

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
    /// Prefetch tiered data with filter settings for a sample and save to file
    Prefetch(PrefetchArgs),
    /// Plot the reference plate layout
    PlotPlate(PlotPlateArgs),
    /// Plot the reference plate review results
    PlotReview(PlotReviewArgs),
    /// Review diagnostic outcome of a review table against the reference data
    Review(ReviewArgs),
    /// Diagnose samples on the reference plate using the generative practitioner
    DiagnoseLocal(DiagnoseLocalArgs),
    /// Debug the pathogen calls made in a set of diagnostic practitioner outputs
    DebugPathogen(DebugPathogenArgs),
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("user")
    .required(false)
    .multiple(true)
    .args(&["controls", "tags"])
)]
#[clap(group = ArgGroup::new("file")
    .args(&["json"])
)]
pub struct EvaluateArgs {
    /// Sample identifier for query
    #[clap(long, short = 's')]
    pub sample: String,
    /// Control sample identifiers associated with sample
    #[clap(long, short = 'c', num_args(0..), group = "user", help_heading = "Validation Query")]
    pub controls: Vec<String>,
    /// Tags for sample and model query
    #[clap(long, short = 't', num_args(0..), group = "user", help_heading = "Validation Query")]
    pub tags: Vec<String>,
    /// JSON file of request schema
    #[clap(long, short = 'j', group = "file", help_heading = "Schema Query")]
    pub json: Option<PathBuf>,
    /// Check for contamination history outliers to be removed from prevalence filter
    #[clap(long)]
    pub contam_history: bool,
}

#[derive(Debug, Args)]
pub struct PlotPlateArgs {
    /// File paths for 
    #[clap(long, short = 's')]
    pub sample: String,
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
    #[clap(long, short = 'o', default_value="review_stats.svg")]
    pub output: PathBuf,
    /// Plot 95% CI interval as grey box behind points (assuming normal distribution)
    #[clap(long, short = 'c')]
    pub ci: bool,
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
    /// The length of the sample to generate (in tokens).
    #[arg(short = 'n', long, default_value_t = 1000)]
    pub sample_len: usize,
    /// The temperature used to generate samples, use 0 for greedy sampling.
    #[arg(long, short='t', default_value_t = 0.8)]
    pub temperature: f64,
    /// GPU device index to run on.
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
    /// Force overwrite output, otherwise skip if exists
    #[clap(long, short = 'f')]
    pub force: bool,
    /// Model directory for generator model(.gguf) and tokenizer file (.json)
    #[clap(long, default_value=".")]
    pub model_dir: PathBuf,
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


#[derive(Debug, Args)]
pub struct PrefetchArgs {
    /// Reference plate file (.json)
    #[clap(long, short = 'p')]
    pub plate: PathBuf,
    /// Output directory for tiered filter data (.prefetch.json) and configuration (.config.json)
    #[clap(long, short = 'o')]
    pub outdir: PathBuf,
    /// Override for primary filter prevalence contamination outliers included (default) or excluded
    #[clap(long, short='c')]
    pub contam_history: Option<bool>,
    /// Force overwrite output, otherwise skip if exists
    #[clap(long, short = 'f')]
    pub force: bool,
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
