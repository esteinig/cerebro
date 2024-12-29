
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
use report::ReportType;


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
    /// Compile a report from configuration file (.json)
    Compile(CompileArgs),
    /// Compile a report from configuration file (.json)
    Template(TemplateArgs),
}
#[cfg(feature = "cli")]
#[derive(Debug, Args)]
pub struct CompileArgs {
    /// Report template
    #[clap(long, short = 'r')]
    pub report: ReportType,
    /// Report template configuration (.json)
    #[clap(long, short = 'c')]
    pub config: PathBuf,
    /// Compile report to PDF (.pdf)
    #[clap(long, short = 'p')]
    pub pdf: Option<PathBuf>,
    /// Compile report to template (.typ)
    #[clap(long, short = 'o')]
    pub template: Option<PathBuf>,
}

#[cfg(feature = "cli")]
#[derive(Debug, Args)]
pub struct TemplateArgs {
    /// Report template
    #[clap(long, short = 'r')]
    pub report: ReportType,
    /// Report configuration template (.json)
    #[clap(long, short = 'c')]
    pub config: PathBuf,
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
            let mut compiler = crate::compiler::CommandLineReportCompiler::new(
                String::from("/"), args.report
            )?;
            
            let report = compiler.report(&args.config)?;

            if let Some(ref path) = args.pdf {
                compiler.pdf(&report, "pdf.typ".to_string(), &path)?;
            }
            if let Some(ref path) = args.template {
                compiler.template(&report, path)?;
            }
        },
        Commands::Template( args ) => {
            match args.report {
                ReportType::PathogenDetection => {
                    crate::report::PathogenDetectionReport::default().to_json_file(&args.config)?;
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