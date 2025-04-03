use std::path::PathBuf;
use clap::{ArgGroup, Args, Parser, Subcommand};

use crate::plate::{DiagnosticOutcome, SampleType};

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
    /// Evaluate a set of samples with a reference template
    Evaluate(EvaluateArgs),
    /// Plot the plate layout and data association
    Plate(PlateArgs),
    /// Review diagnostic outcome of a review table against the reference data
    Review(ReviewArgs),
    /// Diagnose samples on the reference plate using the generative practitioner
    Diagnose(DiagnoseArgs),
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
}

#[derive(Debug, Args)]
pub struct PlateArgs {
}

#[derive(Debug, Args)]
pub struct DiagnoseArgs {
    /// Reference plate file (.json)
    #[clap(long, short = 'r')]
    pub reference: PathBuf,
    /// Output directory for diagnostic results (.json)
    #[clap(long, short = 'o')]
    pub outdir: PathBuf,
    /// OpenAI model for diagnostic queries
    #[clap(long, short = 'm', default_value="gpt-4o-mini")]
    pub model: String,
}


#[derive(Debug, Args)]
pub struct ReviewArgs {
    /// Reference plate file (.json)
    #[clap(long, short = 'r')]
    pub reference: PathBuf,
    /// Plate review file (.tsv)
    #[clap(long, short = 't')]
    pub review: Option<PathBuf>,
    /// Evaluate reference samples by sample type
    #[clap(long, short = 's')]
    pub sample_type: Option<SampleType>,
    /// Set reference results for these samples to 'None'
    #[clap(long, short = 'n', num_args=1..)]
    pub set_none: Option<Vec<String>>,
    /// Handle references with missing orthogonal data and a positive review call
    #[clap(long, short = 't', default_value="indeterminate")]
    pub missing_orthogonal: DiagnosticOutcome,
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
