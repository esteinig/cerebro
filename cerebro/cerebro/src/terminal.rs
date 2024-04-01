use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

/// Cerebro: metagenomic diagnostic for clinical production
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro", version)]
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
    #[clap(subcommand)]
    /// Pipeline processing
    Pipeline(PipelineCommands),
    #[clap(subcommand)]
    /// Report configurations and compiler
    Report(ReportCommands),
    #[clap(subcommand)]
    /// Utility tools for pipeline and stack
    Tools(ToolCommands),
    #[clap(subcommand)]
    /// Stack configuration and deployment
    Stack(StackCommands),
    #[clap(subcommand)]
    /// Application programming interface client
    Api(ApiCommands),
}

/*
==================================
STACK CONFIGURATION AND DEPLOYMENT
==================================
*/

#[derive(Debug, Subcommand)]
pub enum StackCommands {
    /// Deploy a stack configuration
    Deploy(StackDeployArgs)
}   

#[derive(Debug, Args)]
pub struct StackDeployArgs {
    /// Stack configuration file (.toml)
    #[clap(long, short = 'c', env = "CEREBRO_STACK_CONFIG_FILE")]
    pub config: PathBuf,
    /// Deploy for local development with hot-reloads (unsafe in production)
    #[clap(long, short = 'd')]
    pub dev: bool,
    /// Create a `.trigger` file in the deployed repository for manual rebuilds during local development
    /// 
    /// Can be used to modify the trigger file (in the base repository) to manually
    /// launch crate compilation rather than watching for any change in `src` or 
    /// `Cargo.toml` for hot-rebuilds during development
    #[clap(long, short = 't')]
    pub trigger: bool,
    /// Configured stack output directory
    #[clap(long, short = 'o', env = "CEREBRO_STACK_CONFIG_DIR")]
    pub outdir: PathBuf,
    /// Clone the specific branch for this repository
    #[clap(long, short = 'b')]
    pub branch: Option<String>,
    /// Checkout repository at the commit or tag provided, overwrites branch if provided
    #[clap(long, short = 'r')]
    pub revision: Option<String>,
    /// Change the subdomain for multiple concurrent stack deployments 
    /// 
    /// This will configure Traefik to deploy to {app,api}.{subdomain}.{domain}
    /// so that multiple subdomains can be configured on the fly during deployment
    /// for example to `app.dev.cerebro.localhost` or `app.demo.cerebro.localhost`
    #[clap(long, short = 's')]
    pub subdomain: Option<String>,
    /// Public or SSH-like repository URL for cloning into deployment
    #[clap(long, short = 'u', env = "CEREBRO_STACK_GIT_REPO_URL", default_value="git@github.com:esteinig/cerebro.git")]
    pub git_url: String,

    /// Use libgit2 bindings instead of system call (feature not stable)
    #[cfg(feature = "libgit")]
    #[clap(long, short = 'l')]
    pub libgit: bool,
    /// SSH private key for repository clone via SSH
    #[cfg(feature = "libgit")]
    #[clap(long, short = 'k', env = "CEREBRO_STACK_GIT_SSH_PRIVATE_KEY")]
    pub git_ssh_key: Option<PathBuf>,
    /// SSH private key passphrase 
    #[cfg(feature = "libgit")]
    #[clap(long, short = 'p', env = "CEREBRO_STACK_GIT_SSH_PRIVATE_KEY_PWD")]
    pub git_ssh_pwd: Option<String>,
}


/*
===============================
PIPELINE PARSERS AND PROCESSORS
===============================
*/

#[derive(Debug, Subcommand)]
pub enum PipelineCommands {
    /// Parse and process pipeline results
    Process(PipelineProcessArgs),
    /// Quality control table from processed results
    Quality(PipelineQualityArgs),
    /// Create a sample sheet from the input directory
    SampleSheet(PipelineSampleSheetArgs)
}

#[derive(Debug, Args)]
pub struct PipelineProcessArgs {
    /// Pipeline sample results directory for processing
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Taxonomy directory containing 'names.dmp' and 'nodes.dmp' used for classification (NCBI)
    /// 
    /// Must be present for classification processing, otherwise only
    /// the quality control module is proccessed.
    #[clap(long, short = 't')]
    pub taxonomy: Option<PathBuf>,
    /// Output file of processed sample database model (.json)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Optional sample identifier to use for output instead of result directory name
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,

}


#[derive(Debug, Args)]
pub struct PipelineQualityArgs {
    /// Processed pipeline result samples (.json)
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output file for quality control table (.tsv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Include header in quality control table
    #[clap(long, short = 'H')]
    pub header: bool,
    /// Input mass of ERCC or EDCC for biomass calculations (in pg)
    #[clap(long, short = 'e')]
    pub ercc_mass: Option<f64>,

}


#[derive(Debug, Args)]
pub struct PipelineSampleSheetArgs {
    /// Processed pipeline result samples (.json)
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output sample sheet file (.csv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Sample glob - pattern matching to find paired samples
    /// 
    /// The glob string to specify paired sample extensions should be in format: {forward_extension,reverse_extension}
    /// where the wildcard specifies the sample identifier, for example: *_{R1,R2}.fastq.gz will match files
    /// "Sample1_R1.fastq.gz" and "Sample1_R2.fastq.gz" to sample identifier "Sample1"
    #[clap(long, short = 'g', default_value = "*_{R1_001,R2_001}.fastq.gz")]
    pub glob: String,
    /// Run identifier - if not provided uses input directory name
    /// 
    /// If you want to fill in custom run identifiers for each sample, 
    /// you can provide an empty string ("") and edit the sample sheet 
    /// manually.
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Run date - if not provided uses current date (YYYYMMDD)
    /// 
    /// If you want to fill in custom run dates for each sample, 
    /// you can provide an empty string ("") and edit the sample sheet 
    /// manually.
    #[clap(long, short = 'd')]
    pub run_date: Option<String>,
    /// Sample group - if not provided sample group designation is an empty string
    /// 
    /// Sample groups can be specified manually for larger runs containing
    /// sampels from multiple experimental groups - these are later available
    /// in the front-end application
    #[clap(long, short = 's')]
    pub sample_group: Option<String>,
    /// Sample type - if not provided sample type designation is an empty string
    /// 
    /// Sample types can be specified manually for larger runs containing
    /// sampels from multiple biological sources- these are later available
    /// in the front-end application
    #[clap(long, short = 't')]
    pub sample_type: Option<String>,
    /// ERCC input mass in picogram - if not provided input mass is 0
    /// 
    /// In the validation experiments, we test different input masses per sample.
    /// Generally not needed and can be overwritten with options in the fixed
    /// workflow settings later. Set to 25 pg for standard ERCC.
    #[clap(long, short = 'e')]
    pub ercc_input: Option<f64>,
    /// Allow symlink target reading for glob file walking
    #[clap(long, short = 'l')]
    pub symlinks: bool
}

/*
=======================
REPORT TEMPLATE ENGINE
======================
*/

#[derive(Debug, Subcommand)]
pub enum ReportCommands {
    /// Generate clinical reports from data and templates
    Compile(ReportCompileArgs),
}

#[derive(Debug, Args)]
pub struct ReportCompileArgs {
    /// Output report LaTeX template
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Compile LaTeX template into output report
    #[cfg(feature = "pdf")]
    #[clap(long, short = 'p')]
    pub pdf: bool,
    /// Base template configuration file (.toml)
    #[clap(long, short = 'c')]
    pub base_config: PathBuf,
    /// Sample template configuration file (.toml)
    #[clap(long, short = 's')]
    pub sample_config: Option<Vec<PathBuf>>,
    /// Complete: patient template configuration file (.toml)
    #[clap(long, short = 'P')]
    pub patient_config: Option<PathBuf>,
    /// Partial: patient header template configuration file (.toml)
    #[clap(long, short = 'H')]
    pub patient_header_config: Option<PathBuf>,
    /// Partial: patient result template configuration file (.toml)
    #[clap(long, short = 'R')]
    pub patient_result_config: Option<PathBuf>,
}



#[derive(Debug, Subcommand)]
pub enum AnalysisCommands {

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
