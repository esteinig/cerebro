use cerebro_pipe::{modules::quality::{write_quality_tsv, QualityControl}, nextflow::{panviral::PanviralOutput, pathogen::PathogenOutput, quality::{QualityControlFiles, QualityControlOutput}}, terminal::{App, Commands, ProcessCommands, TablesCommands, ToolsCommands}, tools::{scan::ScanReads, sheet::SampleSheet, umi::Umi}, utils::init_logger};
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
                        &args.input, args.id.clone(), args.background
                    )?;

                    QualityControl::from_panviral(&output).to_json(&args.qc)?;
                    
                }

                ProcessCommands::Pathogen(args) => {
                    
                    let output = PathogenOutput::from(
                        &args.input, args.id.clone(), args.background
                    )?;

                    QualityControl::from_pathogen(&output).to_json(&args.qc)?;

                }

                ProcessCommands::Quality(args) => {
                    
                    let output = QualityControlOutput::from(
                        &args.input, args.id.clone(), args.background
                    )?;

                    QualityControl::from_quality(&output).to_json(&args.qc)?;

                }
            }
        },
        Commands::Tools(subcommand) => {
            match subcommand {

                ToolsCommands::ScanReads(args) => {
                    let scanner = ScanReads::new(&args.input);
                    let report = scanner.report()?;
                    
                    match args.json {
                        Some(ref path) => report.to_json(path)?,
                        None => println!("{:#?}", report)
                    }
                },
                ToolsCommands::SampleSheet( args ) => {

                    let sample_sheet = SampleSheet::new(
                        &args.directory, 
                        args.prefix.clone(),
                        &args.glob, 
                        args.ont, 
                        args.run_id.clone(), 
                        args.run_date.clone(),
                        args.sample_group.clone(),
                        args.sample_type.clone(),
                        args.ercc_input,
                        args.symlinks,
                        args.recursive,
                        args.not_unique
                    )?;
                    sample_sheet.write(&args.output)?;
                    
                }
                ToolsCommands::UmiTrim( args ) => {
                    let umi_tool = Umi::new(None, None);
                    umi_tool.trim_umi_index(&args.umi, &args.output, args.trim)?;
                },
                ToolsCommands::UmiCheck( args ) => {
                    let umi_tool = Umi::new(None, None);
                    umi_tool.check_umi_in_read(&args.input, &args.output)?;
                },
                ToolsCommands::UmiPrepCalib( args ) => {
                    let umi_tool = Umi::new(None, None);
                    umi_tool.prepare_calib(&args.input, &args.output)?;
                },
                ToolsCommands::UmiDedupNaive( args ) => {
                    let umi_tool = Umi::new(None, None);
                    let report = umi_tool.naive_dedup(&args.input, &args.output, &args.reads, &args.clusters,  !&args.no_umi, &args.delim, args.head)?;
                    
                    if let Some(path) = &args.json {
                        report.to_json(path)?;
                    }
                },
                ToolsCommands::UmiDedupCalib( args ) => {
                    let umi_tool = Umi::new(None, None);
                    umi_tool.calib_dedup(&args.input, &args.output, &args.clusters, args.summary.clone(), args.identifiers.clone())?;
                },
            }
        },

        Commands::Tables(subcommand) => {
            match subcommand {

                TablesCommands::QualityControl(args) => {
                    write_quality_tsv(&args.json, &args.reads, &args.controls, &args.background)?
                }
            }
        }
        _ => {}
    }


    Ok(())
}