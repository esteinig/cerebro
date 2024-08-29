use cerebro_pipe::{modules::quality::QualityControl, nextflow::panviral::PanviralOutput, terminal::{App, Commands, ProcessCommands}, utils::init_logger};
use clap::Parser;

fn main() -> anyhow::Result<()> {


    init_logger();

    let cli = App::parse();

    match &cli.command {
        
        // Login user for access token
        Commands::Process(subcommand) => {
            match subcommand {

                ProcessCommands::Panviral(args) => {
                    
                    let output = PanviralOutput::from(
                        &args.input, args.id.clone()
                    )?;

                    QualityControl::from_panviral(&output).to_json(&args.qc)?;
                    
                }


            }
        },

        _ => {}
    }


    Ok(())
}