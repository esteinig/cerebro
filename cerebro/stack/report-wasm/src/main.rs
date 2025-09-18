
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
#[cfg(feature = "cli")]
use report::TemplateFormat;
#[cfg(feature = "cli")]
use report::{ReportFormat, ReportType, PathogenDetectionReport, TrainingCompletionReport};
#[cfg(feature = "cli")]
use compiler::LibraryReportCompiler;

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
    /// Report configuration template
    #[clap(long, short = 'c')]
    pub config: PathBuf,
    /// Report configuration template format
    #[clap(long, short = 't', default_value="json")]
    pub template: TemplateFormat,
    /// Output format
    #[clap(long, short = 'f', default_value="pdf")]
    pub format: ReportFormat,
    /// Report type 
    #[clap(long, short = 'r', default_value="pathogen-detection")]
    pub report: ReportType,
    /// Optional logo file (PNG) to replace default logo
    #[clap(long, short = 'l')]
    pub logo: Option<PathBuf>,
    /// Logo width as string compatible with Typst 'png(width: ...)' function otherwise value from configuration template (e.g. 10% or 60mm) 
    #[clap(long, short = 'w')]
    pub width: Option<String>,
}

#[cfg(feature = "cli")]
#[derive(Debug, Args)]
pub struct TemplateArgs {
    /// Report template output (.json | .toml)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Report configuration template format
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

            let mut compiler = LibraryReportCompiler::new(
                String::from("/"), args.report
            )?;
            
            if let Some(logo_path) = args.logo {
                compiler.set_logo_from_path(&logo_path)?;
            }

            if let Some(logo_width) = args.width {
                compiler.set_logo_width(logo_width);
            }

            let report = compiler.report_from_path(&args.config, args.template)?;

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
                },
                ReportType::TrainingCompletion => {
                    match args.template {
                        TemplateFormat::Json => TrainingCompletionReport::default().to_json(&args.output)?,
                        TemplateFormat::Toml => TrainingCompletionReport::default().to_toml(&args.output)?,
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
