use std::path::PathBuf;
use clap::{ArgGroup, Args, Parser, Subcommand};

use crate::gpt::{ClinicalContext, GptModel};


#[cfg(feature = "local")]
use crate::llama::LlamaArgs;

#[cfg(feature = "local")]
use crate::quantized::QuantizedArgs;

#[cfg(feature = "local")]
use crate::qwen::QwenArgs;

#[cfg(feature = "local")]
use crate::text::TextGeneratorArgs;

/// Cerebro: metagenomic generative practitioner (GPT)
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "meta-gp", version)]
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
    /// Diagnose a sample given the evidence and decision tree
    Diagnose(DiagnoseArgs),
    
    #[cfg(feature = "local")]
    /// Run local text generation wrapper
    Generate(TextGeneratorArgs),

    #[cfg(feature = "local")]
    /// Run local text generation with Llama
    Llama(LlamaArgs),

    #[cfg(feature = "local")]
    /// Run local text generation with quantized models
    Quantized(QuantizedArgs),

    #[cfg(feature = "local")]
    /// Run local text generation with quantized models
    Qwen(QwenArgs),
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("user")
    .required(false)
    .multiple(true)
    .args(&["controls", "tags"])
)]
#[clap(group = ArgGroup::new("file")
    .args(&["json"])
    .required(false)
    .multiple(true)
)]
pub struct DiagnoseArgs {
    /// Sample identifier for query
    #[clap(long, short = 's')]
    pub sample: String,
    /// Control sample identifiers associated with sample
    #[clap(long, short = 'c', num_args(0..), group = "user", help_heading = "Validation Query")]
    pub controls: Option<Vec<String>>,
    /// Tags for sample and model query
    #[clap(long, short = 't', num_args(0..), group = "user", help_heading = "Validation Query")]
    pub tags: Option<Vec<String>>,
    /// JSON file of request schema
    #[clap(long, short = 'j', group = "file", help_heading = "Schema Query")]
    pub json: Option<PathBuf>,
    /// OpenAI model for diagnostic queries
    #[clap(long, short = 'm', default_value="3o-mini")]
    pub model: GptModel,
    /// Dry run printing only decision tree and validating inputs
    #[clap(long, short = 'd')]
    pub dry_run: bool,
    /// State log at completion of diagnostic evidence synthesis
    #[clap(long, short = 'l')]
    pub state_log: Option<PathBuf>,
    /// Diagnostic log at completion of diagnostic evidence synthesis
    #[clap(long, short = 'o', default_value="diagnostic.log")]
    pub diagnostic_log: PathBuf,
    /// Clinical context model for diagnostic queries
    #[clap(long, default_value="csf")]
    pub clinical_context: ClinicalContext,
    /// Add memory of the diagostic decision tree to key decision points
    #[clap(long, short = 'd')]
    pub diagnostic_memory: bool,
    /// Check for contamination history outliers to be removed from prevalence filter
    #[clap(long, short = 'h')]
    pub contam_history: bool,
    /// Ignore these taxstr as filter for threshold values
    #[clap(long, short = 'i', num_args=1..)]
    pub ignore_taxstr: Option<Vec<String>>,
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
