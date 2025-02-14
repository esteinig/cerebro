
#[cfg(feature = "cli")]
mod compiler;
#[cfg(feature = "cli")]
mod report;
#[cfg(feature = "cli")]
mod world;

#[cfg(feature = "cli")]
use std::path::PathBuf;
#[cfg(feature = "cli")]
use clap::{Args, Parser, Subcommand};
use report::TemplateFormat;
#[cfg(feature = "cli")]
use report::{ReportFormat, ReportType, PathogenDetectionReport};
#[cfg(feature = "cli")]
use compiler::CommandLineReportCompiler;

#[cfg(feature = "cli")]
/// Cerebro: clinical report compiler for production
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-report", version)]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Commands,
}

#[cfg(feature = "cli")]
#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Compile a report (JSON | TOML)
    Compile(CompileArgs),
    /// Output a report template (JSON | TOML)
    Template(TemplateArgs),
}

#[cfg(feature = "cli")]
#[derive(Debug, Args)]
pub struct CompileArgs {
    /// Report output file
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Report configuration
    #[clap(long, short = 'c')]
    pub config: PathBuf,
    /// Report template format
    #[clap(long, short = 't', default_value="json")]
    pub template: TemplateFormat,
    /// Output format
    #[clap(long, short = 'f', default_value="pdf")]
    pub format: ReportFormat,
    /// Report type 
    #[clap(long, short = 'r', default_value="pathogen-detection")]
    pub report: ReportType,
}

#[cfg(feature = "cli")]
#[derive(Debug, Args)]
pub struct TemplateArgs {
    /// Report template output (.json | .toml)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Report template format
    #[clap(long, short = 't', default_value="json")]
    pub template: TemplateFormat,
    /// Report type
    #[clap(long, short = 'r', default_value="pathogen-detection")]
    pub report: ReportType,
}

#[cfg(feature = "cli")]
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


#[cfg(feature = "cli")]
fn main() -> anyhow::Result<()> {

    let cli = Cli::parse();

    match cli.command {
        Commands::Compile( args ) => {

            let mut compiler = CommandLineReportCompiler::new(
                String::from("/"), args.report
            )?;
            
            let report = compiler.report(&args.config, args.template)?;

            match args.format {
                ReportFormat::Pdf => {
                    compiler.pdf(&report, String::from("/report"), &args.output)?;
                },
                ReportFormat::Svg => {
                    compiler.svg(&report, String::from("/report"), &args.output)?;
                },
                ReportFormat::Typst => {
                    compiler.typst(&report, &args.output)?;
                },
            }
        },
        Commands::Template( args ) => {
            match args.report {
                ReportType::PathogenDetection => {
                    match args.template {
                        TemplateFormat::Json => PathogenDetectionReport::default().to_json(&args.output)?,
                        TemplateFormat::Toml => PathogenDetectionReport::default().to_toml(&args.output)?,
                    }
                }
            }
        }
    }

    Ok(())
}

#[cfg(not(feature = "cli"))]
fn main() {
    eprintln!("Command-line interface feature is not enabled. You can enable it with: `--features cli`.");
}
