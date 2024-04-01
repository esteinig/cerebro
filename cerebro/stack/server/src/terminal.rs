use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

/// Cerebro: production stack server
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-server", version)]
pub struct App {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Run the stack server
    Run(StackServerArgs)
}

#[derive(Debug, Args)]
pub struct StackServerArgs {
    /// Configuration file or environmental variable pointing to configuration file
    #[clap(long, short = 'c', env = "CEREBRO_CONFIG_FILE")]
    pub config: PathBuf,
    /// Host address
    #[clap(long, short = 'H', default_value="127.0.0.1", env = "CEREBRO_SERVER_HOST")]
    pub host: String,
    /// Host port
    #[clap(long, short = 'P', default_value="8080", env = "CEREBRO_SERVER_PORT")]
    pub port: u16,
    /// Threads
    #[clap(long, short = 't', default_value="4")]
    pub threads: usize,
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
