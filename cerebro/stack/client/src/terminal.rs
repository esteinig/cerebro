use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

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
