// Cerebro: metagenomic and -transcriptomic diagnostics for clinical production environments
// Copyright (C) 2023  Eike Steinig

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

use crate::api::report::ClinicalReport;
use crate::stack::watcher;
use clap::Parser;
use rayon::prelude::*;
use pipeline::sheet::SampleSheet;
use stack::watcher::SlackConfig;
use terminal::{ReportCommands, UtilCommands};
use tools::modules::anon::ReadAnonymizer;
use tools::modules::slack::{SlackMessage, SlackMessenger};
use tools::modules::split::Splitter;
use tools::modules::subset::FastaSubset;
use tools::modules::virus::{VirusAlignmentSummary, AnnotationOptions};
use crate::tools::modules::umi::Umi;
use crate::utils::init_logger;

use crate::terminal::{
    App, 
    Commands, 
    PipelineCommands,
    StackCommands, 
    ApiCommands,
    ToolCommands,
    UmiCommands,
    VirusCommands,
    ApiProjectCommands
};

mod api;
mod tools;
mod utils;
mod stack;
mod data;
mod pipeline;
mod terminal;

fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    match &cli.command {
        
        Commands::Pipeline( subcommand ) => {
            match subcommand {
                // Parse and process a pipeline results into a sample model output
                PipelineCommands::Process( args ) => {
                    args.input.par_iter().for_each(|results| {
                        let sample = pipeline::sample::WorkflowSample::new(&results, &args.sample_id, args.taxonomy.clone(), true).expect(
                            &format!("Failed to parse workflow sample directory: {}", results.display())
                        );  // SPECIES-LEVEL
                        sample.write_json(&args.output).expect(
                            &format!("Failed to write sample model to: {}", args.output.display())
                        );
                    });
                },
                // Quality control table
                PipelineCommands::Quality( args ) => {
                    pipeline::utils::create_qc_table(args.input.clone(), &args.output, args.header, args.ercc_mass)?;
                }
                // Sample sheet creation
                PipelineCommands::SampleSheet( args ) => {
                    let sample_sheet = SampleSheet::new(
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
                _ => {}
            }
        },

        Commands::Analysis( subcommand ) => {


        },

        Commands::Report( subcommand ) => {
            match subcommand {
                ReportCommands::Compile( args ) => {
                    let clinical_report = ClinicalReport::from_toml(
                        &args.base_config, 
                        args.patient_config.clone(), 
                        args.patient_header_config.clone(), 
                        args.patient_result_config.clone(), 
                        args.sample_config.clone()
                    )?;

                    #[cfg(feature = "pdf")]
                    {
                        if args.pdf {
                            clinical_report.render_pdf(Some(args.output.clone()), None)?;
                        } else {
                            clinical_report.render_template(Some(args.output.clone()))?;
                        }
                        std::process::exit(0)
                    }
                    clinical_report.render_template(Some(args.output.clone()))?;
                },
                _ => {}
            }
            
        },
        Commands::Tools( subcommand ) => {
            match subcommand {
                ToolCommands::Umi( subcommand ) => {
                    let umi_tool = Umi::new(None, None);
                    match subcommand {
                        UmiCommands::Trim( args ) => {
                            umi_tool.trim_umi_index(&args.umi, &args.output, args.trim)?;
                        },
                        UmiCommands::Check( args ) => {
                            umi_tool.check_umi_in_read(&args.input, &args.output)?;
                        },
                        UmiCommands::PrepareCalib( args ) => {
                            umi_tool.prepare_calib(&args.input, &args.output)?;
                        },
                        UmiCommands::DedupNaive( args ) => {
                            umi_tool.naive_dedup(&args.input, &args.output, &args.reads, &args.summary, !&args.no_umi)?;
                        },
                        UmiCommands::DedupCalib( args ) => {
                            umi_tool.calib_dedup(&args.input, &args.output, &args.clusters, args.summary.clone(), args.identifiers.clone())?;
                        },
                        _ => {},
                    }
                },
                ToolCommands::Virus( subcommand ) => {
                    match subcommand {
                        VirusCommands::ScanRemap( args ) => {
                            let virus_data = VirusAlignmentSummary::new(
                                &args.id,
                                &args.db, 
                                &args.scan, 
                                &args.remap, 
                                &args.consensus, 
                                &args.coverage,
                                &AnnotationOptions::virosaurus()
                            )?;
                            virus_data.write_summary(&args.output, args.header)
                        }
                    }
                },

                ToolCommands::Utils( subcommand ) => {
                    match subcommand {
                        UtilCommands::Subset( args ) => {
                            
                            let fasta_subset = FastaSubset::from_mash(&args.mash, &args.min_identity, &args.min_shared_hashes)?;
                            let seq_ids = fasta_subset.get_record_ids(&args.group_index, &args.group_by, &args.group_sep)?;
                            fasta_subset.subset(seq_ids, &args.fasta, &args.output, &None, &niffler::compression::Level::Six)?;
                        },
                        UtilCommands::Anonymize( args ) => {
                            
                            let anonymizer = ReadAnonymizer::new(None, niffler::compression::Level::Six)?;
                            if args.input.len() == 1 {
                                anonymizer.anonymize_single_end(&args.input[0], &args.output[0])?;
                            } else if args.input.len() == 2 {
                                anonymizer.anonymize_paired_end(&args.input, &args.output, &args.illumina)?;
                            } else {
                                panic!("{}", "Data is not single or paired end")
                            }
                        },
                        UtilCommands::Split( args ) => {
                            
                            let splitter = Splitter::new(&args.outdir, &None, &niffler::compression::Level::Six)?;
                            for fasta in &args.input {
                                splitter.split(&fasta)?
                            }
                        },
                    
                        UtilCommands::Slack( args ) => {
                            let messenger = SlackMessenger::new(&args.token);
                            messenger.send(
                                &SlackMessage::new(&args.channel, &args.message)
                            ).expect("Failed to send message");
                        },
                        UtilCommands::Watcher( args ) => {

                            if let Err(error) = watcher::watch_production(
                                &args.path, 
                                std::time::Duration::from_secs(args.interval),
                                std::time::Duration::from_secs(args.timeout),
                                std::time::Duration::from_secs(args.timeout_interval),
                                SlackConfig { channel: args.slack_channel.clone(), token: args.slack_token.clone() }
                            ) {
                                log::error!("Error: {error:?}");
                            }
                        },
                    }
                }
                _ => {}
            }
            
        },

        // Subcommands for stack configuration and deployment
        Commands::Stack( subcommand ) => {

            match subcommand {
                // Configure the stack deployment by creating the 
                // output directory with templated assets
                StackCommands::Deploy( args ) => {
                    let mut stack = stack::deploy::Stack::from_toml(&args.config)?;
                    stack.configure( &args.outdir, args.dev, args.subdomain.clone(), args.trigger )?;

                    #[cfg(feature = "libgit")]
                    {
                        stack.clone_and_checkout_repository_libgit(
                            &args.git_url, 
                            args.libgit,
                            args.branch.clone(),
                            args.revision.clone(), 
                            args.git_ssh_key.clone(), 
                            args.git_ssh_pwd.clone()
                        )?;
                        std::process::exit(0)  // premature exit ok

                    }
                    stack.clone_and_checkout_repository_process(
                        &args.git_url, 
                        args.branch.clone(),
                        args.revision.clone()
                    )?;
                    
                },
                // Run the stack server with provided configuration
                // Arguments are re-parsed in the asynchroneous routine
                StackCommands::RunServer( _ ) => {
                    api::server::main()?
                },
                
                _ => {}

            }
        },

        Commands::Api( subcommand ) => {
            
            let client = api::client::CerebroClient::new(
                &cli.url,
                &cli.token,
                match subcommand { ApiCommands::Login(_) | ApiCommands::Status(_) => true, _ => false },
                cli.danger_invalid_certificate,
                &cli.token_file
            )?;

            match subcommand {

                // Login user for access token
                ApiCommands::Login( args ) => {
                    client.login_user(&args.email, &args.password)?
                },
                // Ping servers with token
                ApiCommands::Ping(_) => {
                    client.ping_servers()?
                },
                // Ping server without token
                ApiCommands::Status(_) => {
                    client.ping_status()?
                },            
                // Upload sample models to database
                ApiCommands::Upload( args ) => {

                    let cerebro = api::cerebro::model::Cerebro::from_workflow_sample(
                        &pipeline::sample::WorkflowSample::read_json(&args.input)?,
                        &args.sample_sheet,
                        &args.pipeline_config
                    )?;

                    let cerebro = match &args.replace_sample_id {
                        Some(new_sample_id) => cerebro.update_sample_id(
                            &new_sample_id, args.replace_sample_tags.clone()
                        ),
                        None => cerebro
                    };

                    client.upload_model(&cerebro, &args.team_name, &args.project_name, args.db_name.as_ref())?;

                    if let Some(output_file) = &args.output {
                        cerebro.write_json(output_file)?;
                    }
                },
                // Query sample models for taxon summary
                ApiCommands::Taxa( args ) => {

                    client.taxa_summary(
                        &args.team_name, 
                        &args.project_name, 
                        args.db_name.as_ref(),
                        args.filter_config.as_ref(),
                        args.run_ids.clone(),
                        args.sample_ids.clone(),
                        args.workflow_ids.clone(),
                        args.workflow_names.clone(),
                        &args.output
                    )?;

                },
                // Modify (create-delete-update) a project
                ApiCommands::Project( subcommand ) => {
                    match subcommand {

                        // Create new project in a team database
                        ApiProjectCommands::Create( args ) => {
                            client.create_project(
                                &args.team_name, 
                                &args.db_name,
                                &args.project_name, 
                                &args.project_description, 
                            )?
                        },
                        _ => {}
                    }

                },
                _ => {}

            }

        },

    }
    // match &cli.commands {
    //     cli::Commands::RunServer  { .. } => {
    //         api::server::main()?;
    //     }
    //     cli::Commands::Report { pdf, template, base_config, patient_config , patient_header_config, patient_result_config, sample_config} => {
    //         let report = report::ClinicalReport::from_toml(base_config, patient_config.to_owned(),patient_header_config.to_owned(), patient_result_config.to_owned(), sample_config.to_owned() )?;
    //         report.render_pdf(Some(pdf.to_owned()), template.to_owned())?;
    //     }
    //     cli::Commands::HashPassword { password } => {
    //         hash_password(password)
    //     }
    //     cli::Commands::PingApi { } => {
    //         let client = api::client::CerebroClient::new(&cli.url, &cli.token_env, false, cli.danger_accept_invalid_certs, &cli.token_file)?;
    //         client.ping_servers()?;
    //     }
    //     cli::Commands::CreateSampleSheet  { input, output, glob, run_id, run_date, sample_group, symlinks } => {
    //         let sample_sheet = workflow::utils::SampleSheet::new(input, &glob, &false, run_id, run_date, sample_group, symlinks)?;
    //         sample_sheet.write(output)?;
    //     }
    //     cli::Commands::ParseSample  { input, sample_id, taxonomy, output, all_taxa  } => {

    //         let taxonomy = taxonomy::ncbi::load(taxonomy)?;
    //         let wf_sample = workflow::parser::WorkflowSample::new(input, sample_id, &taxonomy, all_taxa)?;

    //         if let Some(output_file) = output { wf_sample.write_json(output_file)?; }

    //     }
    //     cli::Commands::InsertSample  { input, sample_sheet, workflow_config, team_name, project_name, db_name, output } => {

    //         let workflow_sample = workflow::parser::WorkflowSample::read_json(input)?;
    //         let cerebro = api::cerebro::model::Cerebro::from_workflow_sample(&workflow_sample, sample_sheet, workflow_config)?;

    //         let client = api::client::CerebroClient::new(&cli.url, &cli.token_env, false, cli.danger_accept_invalid_certs, &cli.token_file)?;
    //         client.cerebro_insert_model(&cerebro, team_name, project_name, db_name)?;

    //         if let Some(output_file) = output {
    //             cerebro.write_json(output_file)?;
    //         }
    //     }
    //     cli::Commands::Login  { email, password } => {
    //         let client = api::client::CerebroClient::new(&cli.url, &cli.token_env, true, cli.danger_accept_invalid_certs, &cli.token_file)?;
    //         client.login_user(email, password)?;
    //     }
    //     cli::Commands::QcSummary  { input, output, header } => {
    //         workflow::parser::create_qc_table(input, output, header)?;
    //     }
    //     cli::Commands::Contam  { input, output, include_tags, exclude_tags, min_ercc, samples, domains, max_low_complexity , ercc_mass, taxid_mass } => {
    //         let contaminator = analysis::contam::ContaminationDetection::new(input.clone());
    //         let filter = analysis::contam::ContaminationFilter::new(include_tags.clone(), exclude_tags.clone(), *min_ercc, *samples, *max_low_complexity, domains.clone());
    //         contaminator.run(&filter, *ercc_mass, taxid_mass.clone(), output.clone());
    //     }
    //     cli::Commands::GetDbMetadata  { workdir , force } => {
    //         let db_builder = database::DatabaseBuilder::new(workdir, force)?;
    //         db_builder.get_taxonomy_metadata()?;
    //     },
    //     cli::Commands::GetDbData  { workdir , force } => {

    //         let db_builder = database::DatabaseBuilder::new(workdir, force)?;
            
    //         db_builder.get_database_genomes()?;

    //     },
    //     cli::Commands::Subset {
    //         fasta,
    //         mash,
    //         min_identity,
    //         min_shared_hashes,
    //         group_index,
    //         group_by,
    //         group_sep,
    //         output,
    //         output_format,
    //         compression_level,
    //     } => {
    //         match mash {
    //             Some(path) => {
    //                 let fasta_subset = tools::subset::FastaSubset::from_mash(path, min_identity, min_shared_hashes)?;
    //                 let seq_ids = fasta_subset.get_record_ids(group_index, group_by, group_sep)?;
    //                 fasta_subset.subset(seq_ids, fasta, output, output_format, compression_level)?;
    //             },
    //             None => {}
    //         }
    //     },
    //     cli::Commands::Completeness {
    //         input, 
    //         output,
    //         output_format,
    //         compression_level,
    //         min_completeness
    //     } => {
    //         let complete = tools::completeness::Completeness::new(output_format.to_owned(), compression_level.to_owned())?;
    //         complete.completeness(&input, &output, min_completeness)?;
    //     },
    //     cli::Commands::Annotate {
    //         input,
    //         output,
    //         acc2taxid,
    //         ncbi_tax,
    //         exclude
    //     } => {

    //         let taxonomy = match ncbi_tax {
    //             Some(dir) => Some(taxonomy::ncbi::load(dir)?),
    //             None => None
    //         };

    //         let annotate = tools::annotate::Annotate::new(&acc2taxid)?;
    //         annotate.annotate(&input, &output, &exclude, &taxonomy)?;
    //     },
    //     cli::Commands::Anonymize {
    //         input,
    //         output,
    //         output_format,
    //         compression_level,
    //         fake_illumina_header
    //     } => {
            
    //         cli.validate_input_output_combination()?;

    //         let anonymizer = tools::anon::ReadAnonymizer::new(output_format.to_owned(), compression_level.to_owned())?;
    //         if input.len() == 1 {
    //             anonymizer.anonymize_single_end(&input[0], &output[0])?;
    //         } else if input.len() == 2 {
    //             anonymizer.anonymize_paired_end(&input, &output, fake_illumina_header)?;
    //         } else {
    //             panic!("{}", "Data is not single or paired end")
    //         }
                
    //     },
    //     cli::Commands::Split {
    //         input,
    //         outdir,
    //         output_format,
    //         compression_level
    //     } => {
    //         let splitter = tools::split::Splitter::new(outdir, output_format, compression_level)?;
    //         for fasta in input {
    //             splitter.split(&fasta)?
    //         }
    //     }
    // }

    Ok(())
}