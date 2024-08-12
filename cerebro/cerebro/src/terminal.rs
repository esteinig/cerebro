use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

use cerebro_workflow::terminal::Commands as WorkflowCommands;
use cerebro_report::terminal::Commands as ReportCommands;
use cerebro_client::terminal::Commands as ClientCommands;
use cerebro_watcher::terminal::Commands as WatcherCommands;

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
    /// User team name or identifier for requests that require team specification 
    #[clap(
        long, 
        short = 't', 
        env = "CEREBRO_USER_TEAM",
        hide_env_values = true
    )]
    pub team: Option<String>,
    /// SSL certificate verification is ignored [DANGER]
    #[clap(
        long, 
        env = "CEREBRO_DANGER_ACCEPT_INVALID_TLS_CERTIFICATE"
    )]
    /// SeaweedFS master node address
    #[clap(
        long, 
        short = 's', 
        default_value = "http://fs.cerebro.localhost", 
        env = "CEREBRO_FS_URL"
    )]
    pub fs_url: String,
    /// SeaweedFS master node port
    #[clap(
        long, 
        short = 'p',
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
    #[clap(subcommand)]
    /// File system watchers
    Watcher(WatcherCommands),
    #[clap(subcommand)]
    /// Workflow processing
    Workflow(WorkflowCommands),
    #[clap(subcommand)]
    /// API terminal client
    Client(ClientCommands),
    #[clap(subcommand)]
    /// Report compiler
    Report(ReportCommands),
    #[clap(subcommand)]
    /// Stack configuration and deployment
    Stack(StackCommands),
}


/*
==================================
STACK CONFIGURATION AND DEPLOYMENT
==================================
*/

#[derive(Debug, Subcommand)]
pub enum StackCommands {
    /// Deploy a stack configuration
    Deploy(StackDeployArgs),
    /// Hash a password using the stack hash function
    HashPassword(StackHashPasswordArgs)
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
}

#[derive(Debug, Args)]
pub struct StackHashPasswordArgs {
    /// Password to hash
    #[clap(long, short = 'p', env = "CEREBRO_STACK_HASH_PASSWORD")]
    pub password: String
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
