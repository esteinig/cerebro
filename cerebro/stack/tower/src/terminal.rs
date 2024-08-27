use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

/// Cerebro: production file system watcher 
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-watcher", version)]
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
    /// Watch a file path with the provided configuration
    Watch(WatchArgs),
    /// Send a message to the provided Slack channel
    Slack(SlackArgs)
}

#[derive(Debug, Args)]
pub struct WatchArgs {
    /// Registered tower identifier for pulling staged samples
    #[clap(long, short = 'i')]
    pub tower_id: String,

    /// Main script for pipeline
    #[clap(long, short = 'm')]
    pub main: PathBuf,
    /// Main configuration file for pipeline
    #[clap(long, short = 'c')]
    pub config: PathBuf,
    /// Directory for pipeline databases
    #[clap(long, short = 'd')]
    pub databases: PathBuf,

    /// Working directory for pipeline executions
    #[clap(long, short = 'w', default_value=".")]
    pub workdir: PathBuf,
    /// Cleanup pipeline execution directories of successful runs
    #[clap(long)]
    pub cleanup: bool,
}


#[derive(Debug, Args)]
pub struct SlackArgs {
    /// Slack channel name or identifier
    #[clap(long, short = 'c')]
    pub channel: String,
    /// Simple text message to send
    #[clap(long, short = 'm')]
    pub message: String,
    /// Slack API token
    #[clap(long, short = 't', env = "CEREBRO_SLACK_TOKEN", hide_env_values = true)]
    pub token: String,
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
