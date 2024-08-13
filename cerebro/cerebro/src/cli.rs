// Cerebro: metagenomic and -transcriptomic diagnostics for clinical production environments
// Copyright (C) 2024  Eike Steinig

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

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(unreachable_code)]
#![allow(unreachable_patterns)]

use std::fs::create_dir_all;
use std::time::Duration;

use cerebro_client::error::HttpClientError;
use cerebro_model::api::pipelines::schema::RegisterPipelineSchema;
use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;
use cerebro_model::api::watchers::schema::RegisterWatcherSchema;

use cerebro_watcher::utils::WatcherConfigArgs;
use cerebro_watcher::utils::UploadConfigArgs;

use cerebro_workflow::taxon::TaxonThresholdConfig;
use clap::Parser;
use rayon::prelude::*;

use cerebro::tools::password::hash_password;
use cerebro::utils::init_logger;

use cerebro::terminal::{App, Commands, StackCommands};


fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Workflow( subcommand ) => {
            match subcommand {
                cerebro_workflow::terminal::Commands::Tools(subcommand) => {
                    match subcommand {
                        cerebro_workflow::terminal::ToolCommands::ScanRemap( args ) => {
                            let virus_data = cerebro_workflow::virus::VirusAlignmentSummary::new(
                                &args.id,
                                &args.db, 
                                &args.scan, 
                                &args.remap, 
                                &args.consensus, 
                                &args.coverage,
                                &cerebro_workflow::virus::AnnotationOptions::virosaurus()
                            )?;
                            virus_data.write_summary(&args.output, args.header)
                        },
                        cerebro_workflow::terminal::ToolCommands::SubsetFasta( args ) => {
                                            
                            let fasta_subset = cerebro_workflow::internal::subset::FastaSubset::from_mash(&args.mash, &args.min_identity, &args.min_shared_hashes)?;
                            let seq_ids = fasta_subset.get_record_ids(&args.group_index, &args.group_by, &args.group_sep)?;
                            fasta_subset.subset(seq_ids, &args.fasta, &args.output, &None, &niffler::compression::Level::Six)?;
                        },
                        cerebro_workflow::terminal::ToolCommands::Anonymize( args ) => {
                            
                            let anonymizer = cerebro_workflow::internal::anon::ReadAnonymizer::new(None, niffler::compression::Level::Six)?;
                            if args.input.len() == 1 {
                                anonymizer.anonymize_single_end(&args.input[0], &args.output[0])?;
                            } else if args.input.len() == 2 {
                                anonymizer.anonymize_paired_end(&args.input, &args.output, &args.illumina)?;
                            } else {
                                panic!("{}", "Data is not single or paired end")
                            }
                        },
                        cerebro_workflow::terminal::ToolCommands::SplitFasta( args ) => {
                            
                            let splitter = cerebro_workflow::internal::split::Splitter::new(&args.outdir, &None, &niffler::compression::Level::Six)?;
                            for fasta in &args.input {
                                splitter.split(&fasta)?
                            }
                        },
                        cerebro_workflow::terminal::ToolCommands::UmiTrim( args ) => {
                            let umi_tool = cerebro_workflow::internal::umi::Umi::new(None, None);
                            umi_tool.trim_umi_index(&args.umi, &args.output, args.trim)?;
                        },
                        cerebro_workflow::terminal::ToolCommands::UmiCheck( args ) => {
                            let umi_tool = cerebro_workflow::internal::umi::Umi::new(None, None);
                            umi_tool.check_umi_in_read(&args.input, &args.output)?;
                        },
                        cerebro_workflow::terminal::ToolCommands::UmiPrepCalib( args ) => {
                            let umi_tool = cerebro_workflow::internal::umi::Umi::new(None, None);
                            umi_tool.prepare_calib(&args.input, &args.output)?;
                        },
                        cerebro_workflow::terminal::ToolCommands::UmiDedupNaive( args ) => {
                            let umi_tool = cerebro_workflow::internal::umi::Umi::new(None, None);
                            umi_tool.naive_dedup(&args.input, &args.output, &args.reads, &args.summary, !&args.no_umi)?;
                        },
                        cerebro_workflow::terminal::ToolCommands::UmiDedupCalib( args ) => {
                            let umi_tool = cerebro_workflow::internal::umi::Umi::new(None, None);
                            umi_tool.calib_dedup(&args.input, &args.output, &args.clusters, args.summary.clone(), args.identifiers.clone())?;
                        },
                    }
                },
                // Parse and process a pipeline results into a sample model output
                cerebro_workflow::terminal::Commands::Process( args ) => {
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
                cerebro_workflow::terminal::Commands::Quality( args ) => {
                    cerebro_workflow::utils::create_qc_table(args.input.clone(), &args.output, args.header, args.ercc_mass)?;
                },
                // Taxa table
                cerebro_workflow::terminal::Commands::Taxa( args ) => {

                    let threshold_filter = TaxonThresholdConfig::from_args(args);
                    cerebro_workflow::taxon::taxa_summary(args.input.clone(), &args.output, args.sep, args.header, args.extract, args.filter.clone(), &threshold_filter)?;
                },
                // Sample sheet creation
                cerebro_workflow::terminal::Commands::SampleSheet( args ) => {
                    let sample_sheet = cerebro_workflow::sheet::SampleSheet::new(
                        &args.input, 
                        &args.glob, 
                        args.ont, 
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
        },
        Commands::Report( subcommand ) => {
            match subcommand {
                cerebro_report::terminal::Commands::Compile( args ) => {

                    let clinical_report = cerebro_report::report::ClinicalReport::from_toml(
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
        },
        // Subcommands for stack configuration and deployment
        Commands::Stack( subcommand ) => {
            match subcommand {
                StackCommands::Deploy( args ) => {
                    let mut stack = cerebro::stack::deploy::Stack::from_toml(&args.config)?;
                    stack.configure( &args.outdir, args.dev, args.subdomain.clone(), args.trigger )?;

                    stack.clone_and_checkout_repository_process(
                        &args.git_url, 
                        args.branch.clone(),
                        args.revision.clone()
                    )?;
                },
                StackCommands::HashPassword( args ) => {
                    println!("{}", hash_password(&args.password)?);
                },  
            }
        },
        // Subcommands for stack configuration and deployment
        Commands::Watcher( subcommand ) => {
            match subcommand {
                cerebro_watcher::terminal::Commands::Slack( args ) => {
                    let messenger = cerebro_watcher::slack::SlackClient::new(&args.token);
                    messenger.send(
                        &cerebro_watcher::slack::SlackMessage::new(&args.channel, &args.message)
                    ).expect("Failed to send message");
                },
                cerebro_watcher::terminal::Commands::Watch( args ) => {
                            
                    let api_client = cerebro_client::client::CerebroClient::new(
                        &cli.url,
                        cli.token,
                        false,
                        cli.danger_invalid_certificate,
                        cli.token_file,
                        cli.team.clone()
                    )?;
        
                    let fs_client = cerebro_fs::client::FileSystemClient::new(
                        &api_client, 
                        &cli.fs_url, 
                        &cli.fs_port
                    );
        
                    log::info!("Checking status of Cerebro API at {}",  &api_client.url);
                    api_client.ping_servers()?;
            
                    log::info!("Checking status of SeaweedFS master at {}",  &fs_client.fs_url);
                    fs_client.ping_status()?;
        
                    let watcher = cerebro_watcher::watcher::CerebroWatcher::new(
                        cerebro_model::api::watchers::model::ProductionWatcher::from_args(args, &api_client)?, 
                        api_client, 
                        fs_client, 
                        cerebro_fs::client::UploadConfig::from_args(args), 
                        cerebro_watcher::slack::SlackConfig::from_args(args)
                    )?;
        
                    watcher.watch(
                        &args.path, 
                        Duration::from_secs(args.interval),
                        Duration::from_secs(args.timeout),
                        Duration::from_secs(args.timeout_interval),
                    )?;
                    
                },
            }
        },
        Commands::Client( subcommand ) => {


            let client = cerebro_client::client::CerebroClient::new(
                &cli.url,
                cli.token,
                match &subcommand { cerebro_client::terminal::Commands::Login(_) | cerebro_client::terminal::Commands::PingStatus(_) => true, _ => false },
                cli.danger_invalid_certificate,
                cli.token_file,
                cli.team.clone()
            )?;

            match subcommand {
                cerebro_client::terminal::Commands::Login( args ) => {
                    client.login_user(&args.email, args.password.clone())?
                },
                cerebro_client::terminal::Commands::PingServer(_) => {
                    client.ping_servers()?
                },
                cerebro_client::terminal::Commands::PingStatus(_) => {
                    client.ping_status()?
                },     
                cerebro_client::terminal::Commands::Pipeline( subcommand ) => {
                    match subcommand {
                        cerebro_client::terminal::ApiPipelineCommands::Register( args ) => {
                            
                            let register_pipeline_schema = RegisterPipelineSchema::new(
                                &args.name, 
                                &args.location, 
                                args.pipeline.clone()
                            );
        
                            client.register_pipeline(
                                &register_pipeline_schema,
                                true
                            )?;
        
                            if let Some(path) = &args.json {
                                register_pipeline_schema.to_json(path)?;
                            }
                        },
                        cerebro_client::terminal::ApiPipelineCommands::List( args ) => {
                            client.get_pipelines(args.id.clone(), true)?;
                        },
                        cerebro_client::terminal::ApiPipelineCommands::Ping( args ) => {

                            let pipleine_id = match (args.json.clone(), args.id.clone()) {
                                (Some(path), _) => RegisterPipelineSchema::from_json(&path)?.id,
                                (None, Some(id)) => id.to_owned(),
                                (None, None) => return Err(HttpClientError::PipeineIdentifierArgNotFound.into())
                            };

                            client.ping_pipeline(&pipleine_id, true)?;
                        },
                        cerebro_client::terminal::ApiPipelineCommands::Delete( args ) => {
                            
                            let pipleine_id = match (args.json.clone(), args.id.clone()) {
                                (Some(path), _) => RegisterPipelineSchema::from_json(&path)?.id,
                                (None, Some(id)) => id.to_owned(),
                                (None, None) => return Err(HttpClientError::PipeineIdentifierArgNotFound.into())
                            };
        
                            client.delete_pipeline(&pipleine_id)?;
                        },
                    }
                },        
                cerebro_client::terminal::Commands::Watcher(subcommand) => {
                    match subcommand {
                        cerebro_client::terminal::ApiWatcherCommands::Register( args ) => {
                            
                            let register_watcher_schema = RegisterWatcherSchema::new(
                                &args.name, 
                                &args.location,
                                args.format.clone(),
                                args.glob.clone()
                            );
        
                            client.register_watcher(
                                &register_watcher_schema,
                                true
                            )?;
        
                            if let Some(path) = &args.json {
                                register_watcher_schema.to_json(path)?;
                            }
                        },
                        cerebro_client::terminal::ApiWatcherCommands::List( args ) => {
                            client.get_watchers(args.id.clone(), true)?;
                        },
                        cerebro_client::terminal::ApiWatcherCommands::Ping( args ) => {

                            let watcher_id = match (args.json.clone(), args.id.clone()) {
                                (Some(path), _) => RegisterWatcherSchema::from_json(&path)?.id,
                                (None, Some(id)) => id.to_owned(),
                                (None, None) => return Err(HttpClientError::WatcherIdentifierArgNotFound.into())
                            };

                            client.ping_watcher(&watcher_id, true)?;
                        },
                        cerebro_client::terminal::ApiWatcherCommands::Delete( args ) => {
                            
                            let watcher_id = match (args.json.clone(), args.id.clone()) {
                                (Some(path), _) => RegisterWatcherSchema::from_json(&path)?.id,
                                (None, Some(id)) => id.to_owned(),
                                (None, None) => return Err(HttpClientError::WatcherIdentifierArgNotFound.into())
                            };
        
                            client.delete_watcher(&watcher_id)?;
                        },
                    }
                                
                },       
                cerebro_client::terminal::Commands::Stage(subcommand) => {
                    match subcommand {
                        cerebro_client::terminal::ApiStageCommands::Register( args ) => {
                            
                            let schema = RegisterStagedSampleSchema {
                                run_id: args.run_id.clone(),
                                database: args.database.clone(),
                                project: args.project.clone(),
                                pipeline: args.pipeline.clone()
                            };

                            client.register_staged_samples(&schema)?;
                        },
                        cerebro_client::terminal::ApiStageCommands::List( args ) => {
                        },
                        cerebro_client::terminal::ApiStageCommands::Delete( args ) => {
                        },
                    }
                                
                },       
                cerebro_client::terminal::Commands::UploadSample( args ) => {
                    
                    let mut cerebro_models = Vec::new();
                    
                    for path in &args.sample_models {
                        let cerebro = cerebro_model::api::cerebro::model::Cerebro::from_workflow_sample(
                            &cerebro_workflow::sample::WorkflowSample::read_json(&path)?,
                            &args.sample_sheet,
                            &args.pipeline_config
                        )?;

                        let cerebro = match args.replace_sample_id.clone() {
                            Some(new_sample_id) => cerebro.update_sample_id(
                                &new_sample_id, args.replace_sample_tags.clone()
                            ),
                            None => cerebro
                        };

                        cerebro_models.push(cerebro)
                    }
                    
                    for cerebro_model in &cerebro_models {
                        if let Some(model_dir) = &args.model_dir {
                            if !model_dir.exists() || !model_dir.is_dir() {
                                create_dir_all(&model_dir)?
                            }
                            let output_file = model_dir.join(format!("{}.json", cerebro_model.name));
                            log::info!("Writing model to file: {}", output_file.display());
                            cerebro_model.write_json(&output_file)?;
                        }
                    }

                    // client.upload_models(&cerebro_models, &args.team_name, &args.project_name, args.db_name.as_ref())?;
                    
                },
                cerebro_client::terminal::Commands::GetTaxa( args ) => {

                    // client.taxa_summary(
                    //     &args.team_name, 
                    //     &args.project_name, 
                    //     args.db_name.as_ref(),
                    //     args.filter_config.as_ref(),
                    //     args.run_ids.clone(),
                    //     args.sample_ids.clone(),
                    //     args.workflow_ids.clone(),
                    //     args.workflow_names.clone(),
                    //     &args.output
                    // )?;

                },
                cerebro_client::terminal::Commands::GetQuality( args ) => {

                    // client.qc_summary(
                    //     &args.team_name, 
                    //     &args.project_name, 
                    //     args.db_name.as_ref(),
                    //     args.cerebro_ids.clone(),
                    //     args.sample_ids.clone(),
                    //     args.ercc_pg.clone(),
                    //     &args.output
                    // )?;

                },
                cerebro_client::terminal::Commands::Project( subcommand ) => {
                    match subcommand {

                        cerebro_client::terminal::ApiProjectCommands::Create( args ) => {
                            client.create_project(
                                &args.team_name, 
                                &args.db_name,
                                &args.project_name, 
                                &args.project_description, 
                            )?
                        }
                    }
                },
                cerebro_client::terminal::Commands::Team( _subcommand ) => { },
                cerebro_client::terminal::Commands::Database( _subcommand ) => { },
            }
        },
    }
    Ok(())
}