use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

/// Cerebro: production stack server
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-pipeline", version)]
pub struct App {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
      

    #[clap(subcommand)]
    /// Parse and process pipeline results
    Process(ProcessCommands),

    #[clap(subcommand)]
    /// Sample sheet and other input sheet creation
    Sheet(SheetCommands),

    #[clap(subcommand)]
    /// Sample sheet and other input sheet creation
    Table(TableCommands),

    #[clap(subcommand)]
    /// Internal pipeline utilities
    Utils(UtilsCommands),
}




#[derive(Debug, Subcommand)]
pub enum SheetCommands {
}


#[derive(Debug, Subcommand)]
pub enum TableCommands {
}

#[derive(Debug, Subcommand)]
pub enum UtilsCommands {
}


#[derive(Debug, Subcommand)]
pub enum ProcessCommands {
    /// Process panviral enrichment outputs
    Panviral(PanviralArgs),
}


#[derive(Debug, Args)]
pub struct PanviralArgs {
    /// Input directory containing the output files for the panviral pipeline
    #[clap(long, short = 'i', default_value=".")]
    pub input: PathBuf,
    /// Sample identifier to parse, directory basename by default
    #[clap(long, short = 's')]
    pub id: Option<String>,
    /// Output file of processed quality control data
    #[clap(long, short = 'q')]
    pub qc: PathBuf,

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
