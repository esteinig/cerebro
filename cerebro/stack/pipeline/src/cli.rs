use cerebro_pipeline::{modules::{alignment::Alignment, pathogen::{write_pathogen_table, PathogenDetection, PathogenDetectionFilter}, quality::{write_quality_table, QualityControl}}, nextflow::{mag::BlastTaxidMethod, panviral::PanviralOutput, pathogen::PathogenDetectionOutput, quality::{QualityControlFiles, QualityControlOutput}}, terminal::{App, Commands, ProcessCommands, TableCommands, ToolsCommands}, tools::{scan::ScanReads, sheet::SampleSheet, umi::Umi}, utils::init_logger};
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
                        &args.input, 
                        args.input_qc.clone(),
                        args.id.clone(),
                        args.qc_fail_ok
                    )?;
                    let quality_control = QualityControl::from_panviral(&output);
                    let panviral = Alignment::from_panviral(&output)?;

                    if let Some(path) = &args.qc {
                        quality_control.to_json(path)?;
                    }
                    if let Some(path) = &args.panviral {
                        panviral.to_json(path)?;
                    }
                    
                }

                ProcessCommands::Pathogen(args) => {

                    let output = PathogenDetectionOutput::from(
                        &args.input, 
                        args.input_qc.clone(),
                        args.id.clone(),
                        args.taxonomy_directory.clone(),
                        if args.blast_lca { BlastTaxidMethod::LCA } else { BlastTaxidMethod::HighestBitscore },
                        args.qc_fail_ok
                    )?;

                    let quality_control = QualityControl::from_pathogen(&output);

                    let pathogen_detection = PathogenDetection::from_pathogen(
                        &output, 
                        &quality_control, 
                        args.paired_end
                    )?;

                    if let Some(path) = &args.qc {
                        quality_control.to_json(path)?;
                    }
                    if let Some(path) = &args.pathogen {
                        pathogen_detection.to_json(path)?;
                    }
                }

                ProcessCommands::Quality(args) => {
                    let output = QualityControlOutput::from(
                        &args.input, args.id.clone(), args.qc_fail_ok
                    )?;
                    
                    let quality_control = QualityControl::from_quality(&output);
                    
                    if let Some(path) = &args.qc {
                        quality_control.to_json(path)?;
                    }
                }
            }
        },
        Commands::Tools(subcommand) => {
            match subcommand {
                ToolsCommands::Download( args ) => {
                    let dl = args.clone().validate_and_build()?;
        
                    if args.list { dl.list(); } else { dl.download_index()?; }
                }
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
                        args.sample_date.clone(),
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

                    let report = umi_tool.naive_dedup(
                        &args.input, 
                        &args.output, 
                        &args.reads, 
                        &args.clusters,  
                        !&args.no_umi, 
                        &args.delim, 
                        args.head
                    )?;
                    
                    if let Some(path) = &args.json {
                        report.to_json(path)?;
                    }
                },
                ToolsCommands::UmiDedupCalib( args ) => {
                    let umi_tool = Umi::new(None, None);

                    umi_tool.calib_dedup(
                        &args.input, 
                        &args.output, 
                        &args.clusters, 
                        args.summary.clone(), 
                        args.identifiers.clone()
                    )?;
                },
            }
        },

        Commands::Table(subcommand) => {
            match subcommand {

                TableCommands::QualityControl(args) => {
                    write_quality_table(
                        &args.json, 
                        &args.reads, 
                        &args.controls, 
                        &args.background
                    )?
                },

                TableCommands::PathogenDetection(args) => {
                    
                    write_pathogen_table(
                        &args.json,
                        &args.output, 
                        args.taxonomy.clone(),
                    )?
                }

                TableCommands::FilterTable(args) => {
                    
                }
            }
        }
    }


    Ok(())
}