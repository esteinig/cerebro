
use clap::Parser;
use rayon::prelude::*;

use cerebro_workflow::utils;
use cerebro_workflow::terminal;
use cerebro_workflow::internal;
use cerebro_workflow::virus;

use cerebro_workflow::terminal::Commands;
use cerebro_workflow::terminal::ToolCommands;

fn main() -> anyhow::Result<()> {

    utils::init_logger();

    let cli = terminal::App::parse();

    match &cli.command {
        Commands::Tools(subcommand) => {
            match subcommand {
                ToolCommands::ScanRemap( args ) => {
                    let virus_data = virus::VirusAlignmentSummary::new(
                        &args.id,
                        &args.db, 
                        &args.scan, 
                        &args.remap, 
                        &args.consensus, 
                        &args.coverage,
                        &virus::AnnotationOptions::virosaurus()
                    )?;
                    virus_data.write_summary(&args.output, args.header)
                },
                ToolCommands::SubsetFasta( args ) => {
                                    
                    let fasta_subset = internal::subset::FastaSubset::from_mash(&args.mash, &args.min_identity, &args.min_shared_hashes)?;
                    let seq_ids = fasta_subset.get_record_ids(&args.group_index, &args.group_by, &args.group_sep)?;
                    fasta_subset.subset(seq_ids, &args.fasta, &args.output, &None, &niffler::compression::Level::Six)?;
                },
                ToolCommands::Anonymize( args ) => {
                    
                    let anonymizer = internal::anon::ReadAnonymizer::new(None, niffler::compression::Level::Six)?;
                    if args.input.len() == 1 {
                        anonymizer.anonymize_single_end(&args.input[0], &args.output[0])?;
                    } else if args.input.len() == 2 {
                        anonymizer.anonymize_paired_end(&args.input, &args.output, &args.illumina)?;
                    } else {
                        panic!("{}", "Data is not single or paired end")
                    }
                },
                ToolCommands::SplitFasta( args ) => {
                    
                    let splitter = internal::split::Splitter::new(&args.outdir, &None, &niffler::compression::Level::Six)?;
                    for fasta in &args.input {
                        splitter.split(&fasta)?
                    }
                },
                ToolCommands::UmiTrim( args ) => {
                    let umi_tool = internal::umi::Umi::new(None, None);
                    umi_tool.trim_umi_index(&args.umi, &args.output, args.trim)?;
                },
                ToolCommands::UmiCheck( args ) => {
                    let umi_tool = internal::umi::Umi::new(None, None);
                    umi_tool.check_umi_in_read(&args.input, &args.output)?;
                },
                ToolCommands::UmiPrepCalib( args ) => {
                    let umi_tool = internal::umi::Umi::new(None, None);
                    umi_tool.prepare_calib(&args.input, &args.output)?;
                },
                ToolCommands::UmiDedupNaive( args ) => {
                    let umi_tool = internal::umi::Umi::new(None, None);
                    umi_tool.naive_dedup(&args.input, &args.output, &args.reads, &args.summary, !&args.no_umi)?;
                },
                ToolCommands::UmiDedupCalib( args ) => {
                    let umi_tool = internal::umi::Umi::new(None, None);
                    umi_tool.calib_dedup(&args.input, &args.output, &args.clusters, args.summary.clone(), args.identifiers.clone())?;
                },
            }
        },
        // Parse and process a pipeline results into a sample model output
        Commands::Process( args ) => {
            args.input.par_iter().for_each(|results| {
                let sample = cerebro_workflow::sample::WorkflowSample::new(&results, &args.sample_id, args.taxonomy.clone(), true).expect(
                    &format!("Failed to parse workflow sample directory: {}", results.display())
                );  // SPECIES-LEVEL
                sample.write_json(&args.output).expect(
                    &format!("Failed to write sample model to: {}", args.output.display())
                );
            });
        },
        // Quality control table
        Commands::Quality( args ) => {
            cerebro_workflow::utils::create_qc_table(args.input.clone(), &args.output, args.header, args.ercc_mass)?;
        },
        // Sample sheet creation
        Commands::SampleSheet( args ) => {
            let sample_sheet = cerebro_workflow::sheet::SampleSheet::new(
                &args.input, 
                &args.glob, 
                false, 
                args.run_id.clone(), 
                args.run_date.clone(),
                args.sample_group.clone(),
                args.sample_type.clone(),
                args.ercc_input,
                args.symlinks
            )?;
            sample_sheet.write(&args.output)?;
        }
       
    }

    Ok(())

}
   