use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

/// Cerebro: production file system watcher 
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-report", version)]
pub struct App {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Compile a report from configuration file (.toml)
    Compile(CompileArgs),
}

#[derive(Debug, Args)]
pub struct CompileArgs {
    /// Output report Typst template
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Compiled PDF to output instead of template
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
