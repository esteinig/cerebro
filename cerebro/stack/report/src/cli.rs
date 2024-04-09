#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unreachable_code)]



use clap::Parser;

use cerebro_report::utils::init_logger;
use cerebro_report::report::ClinicalReport;
use cerebro_report::terminal::{App, Commands};

// NB: import error arises from including `cli.rs` in `lib.rs` for obvious reasons

fn main() -> anyhow::Result<()> {
    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Compile( args ) => {

            let clinical_report = ClinicalReport::from_toml(
                &args.base_config, 
                args.patient_config.clone(), 
                args.patient_header_config.clone(), 
                args.patient_result_config.clone(), 
                args.sample_config.clone()
            )?;

            if args.pdf {
                clinical_report.render_pdf(Some(args.output.clone()), None, None)?;
            } else {
                clinical_report.render_template(Some(args.output.clone()))?;
            }
            std::process::exit(0)
        }
    }
    Ok(())
}
