use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

/// Cerebro: production file system watcher 
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-watcher", version)]
pub struct App {
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
    /// File path to watch recursively for input folders
    #[clap(long, short = 'p', default_value=".")]
    pub path: PathBuf,
    /// Interval for polling file path recursively in seconds
    #[clap(long, short = 'i', default_value="3")]
    pub interval: u64,
    /// Timeout in seconds to proceed after no further events on input folder
    #[clap(long, short = 't', default_value="10")]
    pub timeout: u64,
    /// Timeout interval for polling input folder recursively in seconds
    #[clap(long, short = 'm', default_value="1")]
    pub timeout_interval: u64,
    /// Slack API token
    #[clap(long, short = 't', env = "CEREBRO_SLACK_TOKEN", hide_env_values = true)]
    pub slack_token: String,
    /// Slack channel
    #[clap(long, short = 'c', env = "CEREBRO_SLACK_CHANNEL", hide_env_values = true)]
    pub slack_channel: String,
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
