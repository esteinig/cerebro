use std::path::PathBuf;
use cerebro_model::api::{watchers::model::WatcherFormat, pipelines::model::Pipeline};
use clap::{ArgGroup, Args, Parser, Subcommand};

/// Cerebro: production stack server
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-client", version)]
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
    /// User team name or identifier for requests that require team specification 
    #[clap(
        long, 
        short = 't', 
        env = "CEREBRO_USER_TEAM",
        hide_env_values = true
    )]
    pub team: Option<String>,
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
    /// Login user for authentication token 
    Login(ApiLoginArgs),
    /// Ping the server as authenticated user 
    PingServer(ApiPingArgs),
    /// Ping the server as unauthenticated user
    PingStatus(ApiStatusArgs),
    /// Upload pipeline outputs to database
    UploadSample(ApiUploadArgs),
    /// Summary of taxa evidence for requested models
    GetTaxa(ApiTaxaArgs),
    /// Summary of quality control for requested models or samples
    GetQuality(ApiQualityArgs),
    
    /// CRUD operations for production watchers
    #[clap(subcommand)]
    Watcher(ApiWatcherCommands),
    /// CRUD operations for production pipelines
    #[clap(subcommand)]
    Pipeline(ApiPipelineCommands),
    /// CRUD operations for staging production files
    #[clap(subcommand)]
    Stage(ApiStageCommands),
    /// CRUD operations for teams
    #[clap(subcommand)]
    Team(ApiTeamCommands),
    #[clap(subcommand)]
    /// CRUD operations for team databases
    Database(ApiDatabaseCommands),
    #[clap(subcommand)]
    /// CRUD operations for team database projects
    Project(ApiProjectCommands),
}

#[derive(Debug, Args)]
pub struct ApiPingArgs {
}

#[derive(Debug, Args)]
pub struct ApiStatusArgs {
}

#[derive(Debug, Args)]
pub struct ApiLoginArgs {
    /// Registered user email
    #[clap(long, short = 'e', env = "CEREBRO_USER_EMAIL")]
    pub email: String,
    /// Registered user password
    #[clap(long, short = 'p', env = "CEREBRO_USER_PASSWORD")]
    pub password: Option<String>,
}


#[derive(Debug, Args)]
pub struct ApiUploadArgs {
    /// Processed pipeline sample models (*.json)
    #[clap(long, short = 'i', num_args(0..))]
    pub sample_models: Vec<PathBuf>,
    /// Pipeline sample sheet (.csv)
    #[clap(long, short = 's')]
    pub sample_sheet: PathBuf,
    /// Pipeline configuration (.json)
    #[clap(long, short = 'c')]
    pub pipeline_config: PathBuf,
    /// Team name for model upload
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Project name for model upload
    #[clap(long, short = 'p')]
    pub project_name: String,
    /// Database name for model upload
    #[clap(long, short = 'd')]
    pub db_name: Option<String>,
    /// Replace sample identifier before upload
    #[clap(long)]
    pub replace_sample_id: Option<String>,
    /// Replace sample tags before upload
    #[clap(long, num_args(0..))]
    pub replace_sample_tags: Option<Vec<String>>,
    /// Output model as file (.json)
    #[clap(long, short = 'o')]
    pub model_dir: Option<PathBuf>,
}


#[derive(Debug, Args)]
pub struct ApiTaxaArgs {
    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Project name for model query
    #[clap(long, short = 'p')]
    pub project_name: String,
    /// Output summary (.csv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: Option<String>,
    /// Filter JSON that satisfies `CerebroFilterConfig` schema
    #[clap(long, short = 'f')]
    pub filter_config: Option<PathBuf>,
    /// Run identifiers to filter results
    #[clap(long, short = 'r', num_args(0..))]
    pub run_ids: Option<Vec<String>>,
    /// Sample identifiers to filter results
    #[clap(long, short = 's', num_args(0..))]
    pub sample_ids: Option<Vec<String>>,
    /// Workflow identifiers to filter results
    #[clap(long, short = 'w', num_args(0..))]
    pub workflow_ids: Option<Vec<String>>,
    /// Workflow mnemnonic names to filter results
    #[clap(long, short = 'n', num_args(0..))]
    pub workflow_names: Option<Vec<String>>,
}


#[derive(Debug, Args)]
pub struct ApiQualityArgs {
    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Project name for model query
    #[clap(long, short = 'p')]
    pub project_name: String,
    /// Output summary (.csv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: Option<String>,
    /// Cerebro database model identifiers to filter results
    #[clap(long, short = 'r', num_args(0..))]
    pub cerebro_ids: Option<Vec<String>>,
    /// Sample identifiers to filter results
    #[clap(long, short = 's', num_args(0..))]
    pub sample_ids: Option<Vec<String>>,
    /// ERCC control input mass for biomass calculations (picogram)
    #[clap(long, short = 'e')]
    pub ercc_pg: Option<f64>,

}

#[derive(Debug, Args)]
pub struct ApiTeamArgs {
    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Project name for model query
    #[clap(long, short = 'p')]
    pub team_descriptions: String,
}

#[derive(Debug, Subcommand)]
pub enum ApiStageCommands {
    /// Register samples for processing with a production pipeline
    Register(ApiStageRegisterArgs),
    /// List registered samples in the staging area
    List(ApiStageListArgs),
    /// Delete registered samples from the staging area
    Delete(ApiStageDeleteArgs),
}


#[derive(Debug, Args)]
pub struct ApiStageRegisterArgs {
    /// Run identifier or name for sample registration
    #[clap(long, short = 'r')]
    pub run_id: String,
    /// Database for pipeline outputs
    #[clap(long, short = 'd')]
    pub database: String,
    /// Database project for pipeline outputs
    #[clap(long, short = 'p')]
    pub project: String,
    /// Cerebro pipeline for registration
    #[clap(long, short = 'i')]
    pub pipeline: Pipeline,
}

#[derive(Debug, Args)]
pub struct ApiStageListArgs {
    /// Staged sample identifier for single record listing
    #[clap(long, short = 'i')]
    pub id: Option<String>,
}

#[derive(Debug, Args)]
pub struct ApiStageDeleteArgs {
    /// Staged sample identifier generated during registration
    #[clap(long, short = 'i', group = "input")]
    pub id: String,
}


#[derive(Debug, Subcommand)]
pub enum ApiPipelineCommands {
    /// Register a new production pipeline
    Register(ApiPipelineRegisterArgs),
    /// List registered production pipelines
    List(ApiPipelineListArgs),
    /// Delete a registered production pipeline
    Delete(ApiPipelineDeleteArgs),
    /// Ping a registered pipeline to update active status
    Ping(ApiPipelinePingArgs),
}

#[derive(Debug, Args)]
pub struct ApiPipelineRegisterArgs {
    /// Cerebro pipeline for registration
    #[clap(long, short = 'p')]
    pub pipeline: Pipeline,
    /// Pipeline name for registration
    #[clap(long, short = 'n')]
    pub name: String,
    /// Pipeline location for registration
    #[clap(long, short = 'l')]
    pub location: String,
    /// Output registration for future reference (.json)
    #[clap(long, short = 'j')]
    pub json: Option<PathBuf>,
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
pub struct ApiPipelinePingArgs {
    /// Pipeline identifier generated during registration
    #[clap(long, short = 'i', group = "input")]
    pub id: Option<String>,
    /// Pipeline registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>
}

#[derive(Debug, Args)]
pub struct ApiPipelineListArgs {
    /// Pipeline identifier generated during registration for single record listing
    #[clap(long, short = 'i')]
    pub id: Option<String>,
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
pub struct ApiPipelineDeleteArgs {
    /// Pipeline identifier generated during registration
    #[clap(long, short = 'i', group = "input")]
    pub id: Option<String>,
    /// Pipeline registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>
}



#[derive(Debug, Subcommand)]
pub enum ApiWatcherCommands {
    /// Register a new production watcher
    Register(ApiWatcherRegisterArgs),
    /// List registered production watchers
    List(ApiWatcherListArgs),
    /// Delete a registered production watcher
    Delete(ApiWatcherDeleteArgs),
    /// Ping a registered watcher to update active status
    Ping(ApiWatcherPingArgs),
}
#[derive(Debug, Args)]
pub struct ApiWatcherRegisterArgs {
    /// Watcher name for registration
    #[clap(long, short = 'n')]
    pub name: String,
    /// Watcher location for registration
    #[clap(long, short = 'l')]
    pub location: String,
    /// Watcher input format 
    #[clap(long, short = 'f')]
    pub format: WatcherFormat,
    /// Fastq file glob when format is 'fastq' 
    #[clap(long, short = 'g')]
    pub glob: Option<String>,
    /// Output registration for future reference (.json)
    #[clap(long, short = 'j')]
    pub json: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct ApiWatcherListArgs {
    /// Watcher identifier generated during registration
    #[clap(long, short = 'i')]
    pub id: Option<String>,
}


#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
pub struct ApiWatcherPingArgs {
    /// Watcher identifier generated during registration for single record listing
    #[clap(long, short = 'i', group = "input")]
    pub id: Option<String>,
    /// Watcher registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
pub struct ApiWatcherDeleteArgs {
    /// Watcher identifier generated during registration
    #[clap(long, short = 'i', group = "input")]
    pub id: Option<String>,
    /// Watcher registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
}


#[derive(Debug, Subcommand)]
pub enum ApiTeamCommands {
    // Create a new team
    Create(ApiTeamCreateArgs)
}


#[derive(Debug, Subcommand)]
pub enum ApiDatabaseCommands {
    // Create a new team
    Create(ApiDatabaseCreateArgs)
}


#[derive(Debug, Subcommand)]
pub enum ApiProjectCommands {
    // Create a new team
    Create(ApiProjectCreateArgs)
}


#[derive(Debug, Args)]
pub struct ApiTeamCreateArgs {
    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Project name for model query
    #[clap(long, short = 'p')]
    pub team_descriptions: String,
}

#[derive(Debug, Args)]
pub struct ApiDatabaseCreateArgs {
    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Project name for model query
    #[clap(long, short = 'p')]
    pub project_name: String,
    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: Option<String>,
}

#[derive(Debug, Args)]
pub struct ApiProjectCreateArgs {
    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: String,
    /// Project name for model query
    #[clap(long, short = 'n')]
    pub project_name: String,
    /// Project name for model query
    #[clap(long, short = 'i')]
    pub project_description: String,
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
