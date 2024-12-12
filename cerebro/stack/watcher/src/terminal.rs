use std::path::PathBuf;
use cerebro_model::api::{towers::model::Pipeline, watchers::model::WatcherFormat};
use clap::{ArgGroup, Args, Parser, Subcommand};

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
#[clap(group = ArgGroup::new("registered")
    .required(false)
    .args(&["id", "json"])
)]
#[clap(group = ArgGroup::new("new")
    .required(false)
    .multiple(true)
    .args(&["name", "location", "format"])
)]
#[clap(group = ArgGroup::new("automatic")
    .required(false)
    .multiple(true)
    .args(&["tower_id", "pipeline"])
)]
pub struct WatchArgs {
    /// File path to watch recursively for input folders
    #[clap(long, short = 'p')]
    pub path: PathBuf,

    /// Watcher identifier generated during registration
    #[clap(long, short = 'i', group = "registered", help_heading = "Registered Watcher")]
    pub id: Option<String>,
    /// Watcher registration record (.json)
    #[clap(long, short = 'j', group = "registered", help_heading = "Registered Watcher")]
    pub json: Option<PathBuf>,

    /// Watcher name to register new watcher
    #[clap(long, short = 'n', group = "new", help_heading = "New Watcher")]
    pub name: Option<String>,
    /// Watcher location to register new watcher
    #[clap(long, short = 'l', group = "new", help_heading = "New Watcher")]
    pub location: Option<String>,
    /// Watcher input format to register new watcher
    #[clap(long, short = 'f', group = "new", help_heading = "New Watcher")]
    pub format: Option<WatcherFormat>,
    /// Optional file glob, overwrites defaults based on format (see long --help)
    /// 
    /// 'fastq' = "*.fastq.gz", 
    /// 'fastq-pe' = "*_{R1,R2}.fastq.gz", 
    /// 'iseq' = "*_{L001_R1_001,L001_R2_001}.fastq.gz", 
    /// 'nextseq' = "*_{R1_001,R2_001}.fastq.gz"
    #[clap(long)]
    pub glob: Option<String>,
    /// Cerebro FS non-default data center for file upload
    #[clap(long)]
    pub data_center: Option<String>,
    /// Cerebro FS non-default replication strategy for file upload
    #[clap(long)]
    pub replication: Option<String>,
    /// Cerebro FS non-default TTL strategy for file upload
    #[clap(long)]
    pub ttl: Option<String>,
    /// Interval for polling file path recursively in seconds
    #[clap(long, default_value="3")]
    pub interval: u64,
    /// Timeout in seconds to proceed after no further events on input folder
    #[clap(long, default_value="10")]
    pub timeout: u64,
    /// Timeout interval for polling input folder recursively in seconds
    #[clap(long, default_value="1")]
    pub timeout_interval: u64,
    /// Slack API token
    #[clap(long, env = "CEREBRO_SLACK_TOKEN", hide_env_values = true)]
    pub slack_token: Option<String>,
    /// Slack channel
    #[clap(long, env = "CEREBRO_SLACK_CHANNEL", hide_env_values = true)]
    pub slack_channel: Option<String>,

    /// Tower identifier for a registered and tower to access staging area
    #[clap(long, group = "automatic", help_heading = "Automatic Launch")]
    pub tower_id: Option<String>,
    /// Pipeline to launch when files are pulled by the tower
    #[clap(long, group = "automatic", help_heading = "Automatic Launch")]
    pub pipeline: Option<Pipeline>,


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
