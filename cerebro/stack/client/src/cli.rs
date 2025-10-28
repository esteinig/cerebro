use std::fs::create_dir_all;

use cerebro_client::error::HttpClientError;
use cerebro_model::api::cerebro::schema::{PathogenDetectionTableSchema, QualityControlTableSchema};
use cerebro_model::api::jobs::schema::{EnqueueJobRequest, ScheduleJobRequest};
use cerebro_model::api::stage::model::StagedSample;
use cerebro_model::api::towers::schema::RegisterTowerSchema;
use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;
use cerebro_model::api::watchers::schema::RegisterWatcherSchema;
use cerebro_pipeline::modules::alignment::Alignment;
use cerebro_pipeline::modules::pathogen::PathogenDetection;
use cerebro_pipeline::modules::quality::QualityControl;
use clap::Parser;

use cerebro_client::utils::init_logger;
use cerebro_client::client::CerebroClient;
use cerebro_client::terminal::{App, Commands, DatabaseCommands, JobsCommands, ProjectCommands, StageCommands, TowerCommands, TrainingCommands, WatcherCommands};

use cerebro_model::api::cerebro::model::Cerebro;
use cerebro_pipeline::taxa::taxon::TaxonExtraction;


fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    let client = CerebroClient::new(
        &cli.url,
        cli.token,
        match &cli.command { Commands::Login(_) | Commands::PingStatus(_) => true, _ => false },
        cli.danger_invalid_certificate,
        cli.token_file,
        cli.team.clone(),
        cli.db.clone(),
        cli.project.clone()
    )?;

    match &cli.command {
        
        // Login user for access token
        Commands::Login( args ) => {
            client.login_user(&args.email, args.password.clone())?
        },
        // Ping servers with token
        Commands::PingServer(_) => {
            client.ping_servers()?
        },
        // Ping server without token
        Commands::PingStatus(_) => {
            client.ping_status()?
        },            
        // Process and upload sample model to database
        Commands::UploadPathogen( args ) => {
            
            let quality = QualityControl::from_json(&args.quality)?;
            let pathogen = PathogenDetection::from_json(&args.pathogen)?;

            if quality.id != pathogen.id {
                return Err(HttpClientError::PathogenIdentifiersNotMatched(quality.id, pathogen.id).into())
            }

            
            let (team, database, project, run_id) = match &args.stage_json {
                Some(path) => {
                    let staged = StagedSample::from_json(path)?;
                    (staged.team, staged.database, staged.project, staged.run_id)
                },
                None => {
                    match (cli.team.clone(), cli.db.clone(), cli.project.clone(), args.run_id.clone()) {
                        (Some(team), Some(database), Some(project), Some(run_id)) => (team, database, project, run_id),
                        _ => return Err(HttpClientError::MissingUploadParameters.into())
                    }
                }
            };

            let cerebro = Cerebro::from(
                &quality.id,
                &quality,
                &pathogen.get_taxa(
                    &args.taxonomy, 
                    args.strict,
                    args.gtdb_break_monophyly
                )?,
                args.sample_sheet.clone(),
                args.pipeline_config.clone(),
                Some(run_id),
                args.run_date.clone() // TODO: add run configuration to staged sample
            )?;

            
            if let Some(model) = &args.model {
                cerebro.write_json(&model)?;
            }


            let mut client = client.clone();
            client.set_data_auth(&team, &database, &project);  // can be from staged sample

            client.upload_models(
                &Vec::from([cerebro])
            )?;
            
        },
        // Process and upload sample model to database
        Commands::CreatePathogen( args ) => {

            let quality = QualityControl::from_json(&args.quality)?;
            let pathogen = PathogenDetection::from_json(&args.pathogen)?;

            if quality.id != pathogen.id {
                return Err(HttpClientError::PathogenIdentifiersNotMatched(quality.id, pathogen.id).into())
            }

            let cerebro = Cerebro::from(
                &quality.id,
                &quality,
                &pathogen.get_taxa(
                    &args.taxonomy, 
                    args.strict,
                    args.gtdb_break_monophyly
                )?,
                args.sample_sheet.clone(),
                args.pipeline_config.clone(),
                args.run_id.clone(),
                args.run_date.clone()
            )?;

            cerebro.write_json(&args.model)?;

        }
        // Process and upload sample model to database
        Commands::UploadPanviral( args ) => {
            
            let quality = QualityControl::from_json(&args.quality)?;
            let panviral = Alignment::from_json(&args.panviral)?;

            if quality.id != panviral.id {
                return Err(HttpClientError::PanviralIdentifiersNotMatched(quality.id, panviral.id).into())
            }

            let (team, database, project, run_id) = match &args.stage_json {
                Some(path) => {
                    let staged = StagedSample::from_json(path)?;
                    (staged.team, staged.database, staged.project, staged.run_id)
                },
                None => {
                    match (cli.team.clone(), cli.db.clone(), cli.project.clone(), args.run_id.clone()) {
                        (Some(team), Some(database), Some(project), Some(run_id)) => (team, database, project, run_id),
                        _ => return Err(HttpClientError::MissingUploadParameters.into())
                    }
                }
            };

            let cerebro = Cerebro::from(
                &quality.id,
                &quality,
                &panviral.get_taxa(
                    &args.taxonomy, 
                    args.strict,
                    args.gtdb_break_monophyly
                )?,
                args.sample_sheet.clone(),
                args.pipeline_config.clone(),
                Some(run_id),
                args.run_date.clone()
            )?;

            if let Some(model) = &args.model {
                cerebro.write_json(&model)?;
            }
                
            let mut client = client.clone();
            client.set_data_auth(&team, &database, &project); // can be from staged sample

            client.upload_models(
                &Vec::from([cerebro])
            )?;
            
        },
        // Process and upload sample model to database
        Commands::CreatePanviral( args ) => {

            let quality = QualityControl::from_json(&args.quality)?;
            let panviral = Alignment::from_json(&args.panviral)?;

            if quality.id != panviral.id {
                return Err(HttpClientError::PanviralIdentifiersNotMatched(quality.id, panviral.id).into())
            }

            let cerebro = Cerebro::from(
                &quality.id,
                &quality,
                &panviral.get_taxa(
                    &args.taxonomy, 
                    args.strict,
                    args.gtdb_break_monophyly
                )?,
                args.sample_sheet.clone(),
                args.pipeline_config.clone(),
                args.run_id.clone(),
                args.run_date.clone()
            )?;

            cerebro.write_json(&args.model)?;

        }
        Commands::UploadModels( args ) => {
            for path in &args.models {
                log::info!("Reading model from: {}", path.display());
                let model = Cerebro::from_json(path)?;
                client.upload_models(&vec![model])?;
            }
        }
        Commands::DownloadModels( args ) => {

            create_dir_all(&args.outdir)?;
            client.download_models(&args.outdir)?;
        }

        Commands::UpdateModels( args ) => {
            
            client.update_models(
                args.run_tsv.clone()
            )?;
        }

        Commands::DeleteModels( args ) => {
            
            client.delete_models(
                args.name_tsv.clone()
            )?;
        }
        Commands::Tower(subcommand) => {
            match subcommand {
                TowerCommands::Register( args ) => {
                    
                    let schema = RegisterTowerSchema::new(
                        &args.name, 
                        &args.location, 
                        args.pipelines.clone()
                    );

                    client.register_tower(
                        &schema,
                        true
                    )?;

                    if let Some(path) = &args.json {
                        schema.to_json(path)?;
                    }
                },
                TowerCommands::List( args ) => {
                    client.list_towers(args.id.clone(), true)?;
                },
                TowerCommands::Ping( args ) => {

                    let tower_id = match (args.json.clone(), args.id.clone()) {
                        (Some(path), _) => RegisterTowerSchema::from_json(&path)?.id,
                        (None, Some(id)) => id.to_owned(),
                        _ => return Err(HttpClientError::TowerIdentifierNotFound.into())
                    };

                    client.ping_tower(&tower_id,  true)?;
                },
                TowerCommands::Delete( args ) => {
                    if args.all {
                        let confirmation = dialoguer::Confirm::new()
                            .with_prompt("Do you want to delete ALL tower registrations?")
                            .interact()
                            .unwrap();

                        if confirmation {
                            client.delete_tower(
                                None, 
                                None, 
                                None, 
                                None
                            )?;
                        }
                    } else {
                        client.delete_tower(
                            args.id.clone(), 
                            args.json.clone(), 
                            args.name.clone(), 
                            args.location.clone()
                        )?;
                    }
                },
            }
                        
        },
        Commands::Watcher(subcommand) => {
            match subcommand {
                WatcherCommands::Register( args ) => {
                    
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
                WatcherCommands::List( args ) => {
                    client.list_watchers(args.id.clone(), true)?;
                },

                WatcherCommands::Ping( args ) => {

                    let watcher_id = match (args.json.clone(), args.id.clone()) {
                        (Some(path), _) => RegisterWatcherSchema::from_json(&path)?.id,
                        (None, Some(id)) => id.to_owned(),
                        (None, None) => return Err(HttpClientError::WatcherIdentifierArgNotFound.into())
                    };

                    client.ping_watcher(&watcher_id, true)?;
                },
                WatcherCommands::Delete( args ) => {
                    
                    if args.all {
                        let confirmation = dialoguer::Confirm::new()
                            .with_prompt("Do you want to delete ALL watchers?")
                            .interact()
                            .unwrap();

                        if confirmation {
                            client.delete_watcher(
                                None, 
                                None, 
                                None, 
                                None
                            )?;
                        }
                    } else {
                        client.delete_watcher(
                            args.id.clone(), 
                            args.json.clone(), 
                            args.name.clone(), 
                            args.location.clone()
                        )?;
                    }

                },
            }
                        
        },
        Commands::Stage(subcommand) => {
            match subcommand {
                StageCommands::Register( args ) => {
                    
                    let schema = match (args.json.clone(), args.tower_id.clone()) {
                        (Some(path), _) => RegisterStagedSampleSchema::from_tower_json(
                            &path, args.pipeline.clone(), args.file_id.clone(), args.run_id.clone()
                        )?,
                        (None, Some(id)) => RegisterStagedSampleSchema::new(
                            &id, args.pipeline.clone(), args.file_id.clone(), args.run_id.clone()
                        ),
                        (None, None) => return Err(HttpClientError::TowerIdentifierNotFound.into())
                    };

                    client.register_staged_samples(&schema)?;

                },
                StageCommands::List( args ) => {


                    let tower_id = match (args.json.clone(), args.tower_id.clone()) {
                        (Some(path), _) => RegisterTowerSchema::from_json(&path)?.id,
                        (None, Some(id)) => id,
                        (None, None) => return Err(HttpClientError::TowerIdentifierNotFound.into())
                    };

                    client.list_staged_samples(
                        &tower_id, 
                        args.run_id.clone(), 
                        args.sample_id.clone(),
                        true
                    )?;
                },
                StageCommands::Pull( args ) => {

                    let tower_id = match (args.json.clone(), args.tower_id.clone()) {
                        (Some(path), _) => RegisterTowerSchema::from_json(&path)?.id,
                        (None, Some(id)) => id,
                        (None, None) => return Err(HttpClientError::TowerIdentifierNotFound.into())
                    };

                    client.pull_staged_samples(
                        &tower_id, 
                        args.run_id.clone(), 
                        args.sample_id.clone(),
                        &args.outdir,
                        args.delete,
                    )?;

                },
                StageCommands::Delete( args ) => {

                    let tower_id = match (args.json.clone(), args.tower_id.clone()) {
                        (Some(path), _) => RegisterTowerSchema::from_json(&path)?.id,
                        (None, Some(id)) => id,
                        (None, None) => return Err(HttpClientError::TowerIdentifierNotFound.into())
                    };

                    if args.all {
                        let confirmation = dialoguer::Confirm::new()
                            .with_prompt("Do you want to delete ALL staged samples for this staging area?")
                            .interact()
                            .unwrap();

                        if confirmation {
                            client.delete_staged_sample(
                                &tower_id, 
                                None, 
                                None, 
                                None
                            )?;
                        }
                    } else {
                        client.delete_staged_sample(
                            &tower_id, 
                            args.staged_id.clone(), 
                            args.run_id.clone(), 
                            args.sample_id.clone()
                        )?;
                    }                    
                },
            }
                        
        },
        // Query sample models for taxon history and taxon vs host regression analysis
        Commands::GetTaxonHistory( args ) => {

            client.get_taxon_history(
                args.taxon_label.clone(), 
                args.host_label.clone(), 
                args.regression, 
                true,
                args.plot.clone()
            )?;

        },
        // Query sample summaries for quality control data and run/sample/workflow configs
        Commands::GetQualityControlTable( args ) => {

            let schema = QualityControlTableSchema::new(
                args.sample_ids.clone().unwrap_or_default()
            );

            let sample_summaries = client.get_quality_control(
                &schema,
                args.output.clone()
            )?;

            for sample_summary in sample_summaries {
                log::info!("{sample_summary:#?}")
            }

        },

        // Query sample summaries for quality control data and run/sample/workflow configs
        Commands::GetPathogenDetectionTable( args ) => {

            let schema = PathogenDetectionTableSchema::new(
                args.sample_ids.clone().unwrap_or_default(), args.collapse_variants
            );

            client.get_pathogen_detection(
                &schema,
                args.output.clone()
            )?;

        },
        // Training collections for prefetch data
        Commands::Training( subcommand ) => {
            match subcommand {

                // Upload a prefetch entry into the team training database
                TrainingCommands::Upload( args ) => {
                    client.upload_training_prefetch(&args.prefetch, &args.collection, &args.description)?;
                },
                // List prefetch data in a team training database
                TrainingCommands::ListDatasets( args ) => {
                    client.list_training_collections(args.collection.clone())?;
                },
                // Delete prefetch data from the team training database
                TrainingCommands::DeleteDataset( args ) => {
                    
                    if args.name.is_none() {
                        let confirmation = dialoguer::Confirm::new()
                            .with_prompt("Do you want to delete ALL entries for this training collection?")
                            .interact()
                            .unwrap();

                        if confirmation {
                            client.delete_training_data(
                                args.collection.clone(), 
                                args.name.clone()
                            )?;
                        }
                    } else {
                        client.delete_training_data(
                            args.collection.clone(), 
                            args.name.clone()
                        )?;
                    }
                },

                // List session data in a team training database
                TrainingCommands::ListSessions( args ) => {
                    client.list_training_sessions(args.completed, args.user_name.clone(), args.outdir.clone())?;
                },
                // Delete session data in a team training database
                TrainingCommands::DeleteSession( args ) => {
                    client.delete_training_sessions(args.session_id.clone(), args.user_id.clone(), args.incomplete)?;
                },
            }
        },
        // Modify (create-delete-update) a project
        Commands::Project( subcommand ) => {
            match subcommand {

                // Create new project in a team database
                ProjectCommands::Create( args ) => {
                    client.create_project(
                        &args.name, 
                        &args.description, 
                    )?
                }
            }
        },
        Commands::Team( _subcommand ) => { },
        Commands::Database( subcommand ) => {
            match subcommand {

                // Create new project in a team database
                DatabaseCommands::Create( args ) => {
                    client.create_database(
                        &args.name, 
                        &args.description, 
                    )?
                }
            }
        },
        Commands::Jobs(sub) => {
            match sub {
                JobsCommands::LaunchJob { kind, args, queue, retry, reserve_for_seconds } => {
                    
                    // parse JSON string to Value
                    let args_json: serde_json::Value = serde_json::from_str(args)
                        .map_err(|e| HttpClientError::DataResponseFailure(
                            reqwest::StatusCode::BAD_REQUEST, format!("Invalid --args JSON: {e}")
                        ))?;

                    let req = EnqueueJobRequest {
                        kind: kind.clone(),
                        args: args_json,
                        queue: queue.clone(),
                        retry: *retry,
                        reserve_for_seconds: *reserve_for_seconds,
                    };
                    let id = client.enqueue_job(&req)?;
                    println!("{}", id.unwrap_or("error".to_string()));
                }
                JobsCommands::ScheduleJob {
                    kind, 
                    args, 
                    queue, 
                    run_at, 
                    interval_seconds, 
                    retry, 
                    reserve_for_seconds, 
                    enabled
                } => {

                    let args_json: serde_json::Value = serde_json::from_str(args)
                        .map_err(|e| HttpClientError::DataResponseFailure(
                            reqwest::StatusCode::BAD_REQUEST, format!("Invalid --args JSON: {e}")
                        ))?;

                    let run_at_dt = chrono::DateTime::parse_from_rfc3339(run_at)
                        .map_err(|e| HttpClientError::DataResponseFailure(
                            reqwest::StatusCode::BAD_REQUEST, format!("Invalid --run-at RFC3339: {e}")
                        ))?
                        .with_timezone(&chrono::Utc);

                    let req = ScheduleJobRequest {
                        kind: kind.clone(),
                        args: args_json,
                        queue: queue.clone(),
                        run_at: run_at_dt,
                        interval_seconds: *interval_seconds,
                        retry: *retry,
                        reserve_for_seconds: *reserve_for_seconds,
                        enabled: *enabled,
                    };
                    
                    let schedule_id = client.schedule_job(&req)?;
                    println!("{schedule_id}");
                }

                JobsCommands::ScheduleSummary => {
                    let s = client.schedule_status()?;
                    println!("{s}");
                }
                JobsCommands::JobStatus { id } => {
                    let s = client.job_status(id)?;
                    println!("{s}");
                }
                JobsCommands::JobComplete { id } => {
                    let yn = client.job_complete(id)?;
                    println!("{yn}");
                }
            }
        }
    }

    Ok(())

}
   