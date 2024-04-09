
use clap::Parser;

use cerebro_client::utils::init_logger;
use cerebro_client::client::CerebroClient;
use cerebro_client::terminal::{App, Commands, ApiProjectCommands};

use cerebro_workflow::sample::WorkflowSample;
use cerebro_model::api::cerebro::model::Cerebro;

fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    let client = CerebroClient::new(
        &cli.url,
        &cli.token,
        match &cli.command { Commands::Login(_) | Commands::PingStatus(_) => true, _ => false },
        cli.danger_invalid_certificate,
        &cli.token_file
    )?;

    match &cli.command {
        


        // Login user for access token
        Commands::Login( args ) => {
            client.login_user(&args.email, &args.password)?
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
            

            client.upload_models(&cerebro_models, &args.team_name, &args.project_name, args.db_name.as_ref())?;
            
            for cerebro_model in cerebro_models {
                if let Some(output_file) = &args.model_output {
                    log::info!("Writing model to file: {}", output_file.display());
                    cerebro_model.write_json(output_file)?;
                }
            }
            
        },
        // Query sample models for taxon summary
        Commands::GetTaxa( args ) => {

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
        // Query sample models for quality control summary
        Commands::GetQuality( args ) => {

            client.qc_summary(
                &args.team_name, 
                &args.project_name, 
                args.db_name.as_ref(),
                args.cerebro_ids.clone(),
                args.sample_ids.clone(),
                args.ercc_pg.clone(),
                &args.output
            )?;

        },
        // Modify (create-delete-update) a project
        Commands::Project( subcommand ) => {
            match subcommand {

                // Create new project in a team database
                ApiProjectCommands::Create( args ) => {
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
   