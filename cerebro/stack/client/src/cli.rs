
use std::fs::create_dir_all;

use cerebro_client::error::HttpClientError;
use cerebro_model::api::pipelines::schema::RegisterPipelineSchema;
use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;
use cerebro_model::api::watchers::schema::RegisterWatcherSchema;
use clap::Parser;

use cerebro_client::utils::init_logger;
use cerebro_client::client::CerebroClient;
use cerebro_client::terminal::{PipelineCommands, ProjectCommands, StageCommands, WatcherCommands, App, Commands};

use cerebro_workflow::sample::WorkflowSample;
use cerebro_model::api::cerebro::model::Cerebro;

fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    let client = CerebroClient::new(
        &cli.url,
        cli.token,
        match &cli.command { Commands::Login(_) | Commands::PingStatus(_) => true, _ => false },
        cli.danger_invalid_certificate,
        cli.token_file,
        cli.team,
        cli.db,
        cli.project
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
        // Upload sample models to database
        Commands::UploadSample( args ) => {
            
            let mut cerebro_models = Vec::new();
            
            for path in &args.sample_models {
                let cerebro = Cerebro::from_workflow_sample(
                    &WorkflowSample::read_json(&path)?,
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
        Commands::Pipeline(subcommand) => {
            match subcommand {
                PipelineCommands::Register( args ) => {
                    
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
                PipelineCommands::List( args ) => {
                    client.list_pipelines(args.id.clone(), true)?;
                },
                PipelineCommands::Ping( args ) => {

                    let pipleine_id = match (args.json.clone(), args.id.clone()) {
                        (Some(path), _) => RegisterPipelineSchema::from_json(&path)?.id,
                        (None, Some(id)) => id.to_owned(),
                        _ => return Err(HttpClientError::PipelineIdentifierArgNotFound.into())
                    };

                    client.ping_pipeline(&pipleine_id,  true)?;
                },
                PipelineCommands::Delete( args ) => {
                    if args.all {
                        let confirmation = dialoguer::Confirm::new()
                            .with_prompt("Do you want to delete ALL pipelines?")
                            .interact()
                            .unwrap();

                        if confirmation {
                            client.delete_pipeline(
                                None, 
                                None, 
                                None, 
                                None
                            )?;
                        }
                    } else {
                        client.delete_pipeline(
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
                    
                    let schema = match (args.json.clone(), args.id.clone()) {
                        (Some(path), _) => RegisterStagedSampleSchema::from_pipeline_json(
                            &path, args.file_id.clone(), args.run_id.clone()
                        )?,
                        (None, Some(id)) => RegisterStagedSampleSchema::new(
                            &id, args.file_id.clone(), args.run_id.clone()
                        ),
                        (None, None) => return Err(HttpClientError::PipelineIdentifierArgNotFound.into())
                    };

                    client.register_staged_samples(&schema)?;

                },
                StageCommands::List( args ) => {


                    let pipeline_id = match (args.json.clone(), args.id.clone()) {
                        (Some(path), _) => RegisterPipelineSchema::from_json(&path)?.id,
                        (None, Some(id)) => id,
                        (None, None) => return Err(HttpClientError::PipelineIdentifierArgNotFound.into())
                    };

                    client.list_staged_samples(
                        pipeline_id, 
                        args.run_id.clone(), 
                        args.sample_id.clone(),
                        true
                    )?;
                },
                StageCommands::Pull( args ) => {

                    let pipeline_id = match (args.json.clone(), args.id.clone()) {
                        (Some(path), _) => RegisterPipelineSchema::from_json(&path)?.id,
                        (None, Some(id)) => id,
                        (None, None) => return Err(HttpClientError::PipelineIdentifierArgNotFound.into())
                    };


                },
                StageCommands::Delete( args ) => {

                    let pipeline_id = match (args.json.clone(), args.id.clone()) {
                        (Some(path), _) => RegisterPipelineSchema::from_json(&path)?.id,
                        (None, Some(id)) => id,
                        (None, None) => return Err(HttpClientError::PipelineIdentifierArgNotFound.into())
                    };

                    if args.all {
                        let confirmation = dialoguer::Confirm::new()
                            .with_prompt("Do you want to delete ALL staged samples for this staging area?")
                            .interact()
                            .unwrap();

                        if confirmation {
                            client.delete_staged_sample(
                                &pipeline_id, 
                                None, 
                                None, 
                                None
                            )?;
                        }
                    } else {
                        client.delete_staged_sample(
                            &pipeline_id, 
                            args.staged_id.clone(), 
                            args.run_id.clone(), 
                            args.sample_id.clone()
                        )?;
                    }                    
                },
            }
                        
        },
        // Query sample models for taxon summary
        Commands::GetTaxa( args ) => {

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
        // Query sample models for quality control summary
        Commands::GetQuality( args ) => {

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
        // Modify (create-delete-update) a project
        Commands::Project( subcommand ) => {
            match subcommand {

                // Create new project in a team database
                ProjectCommands::Create( args ) => {
                    client.create_project(
                        &args.team_name, 
                        &args.db_name,
                        &args.project_name, 
                        &args.project_description, 
                    )?
                }
            }
        },
        Commands::Team( _subcommand ) => { },
        Commands::Database( _subcommand ) => { },
    }

    Ok(())

}
   