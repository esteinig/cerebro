
use cerebro_model::api::cerebro::schema::GpConfig;
use cerebro_pipeline::taxa::filter::PrevalenceContaminationConfig;
use clap::Parser;
use cerebro_gp::{gpt::{DataBackground, DiagnosticAgent}, terminal::{App, Commands}, utils::init_logger};
use cerebro_client::client::CerebroClient;

#[tokio::main]
async fn main() -> anyhow::Result<(), anyhow::Error> {
    
    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Diagnose( args ) => {

            let api_client = CerebroClient::new(
                &cli.url,
                cli.token,
                false,
                cli.danger_invalid_certificate,
                cli.token_file,
                cli.team,
                cli.db,
                cli.project
            )?;

            log::info!("Checking status of Cerebro API at {}",  &api_client.url);
            api_client.ping_servers()?;
                        
            let mut agent = DiagnosticAgent::new(
                api_client, 
                args.model.clone()
            ).await?;

            // Optionally, print the knowledge graph
            agent.print_decision_tree(&agent.tree, "start");

            if !args.dry_run {

                // Parse generative practicioner configuration
                let mut gp_config = match &args.json {
                    Some(path) => GpConfig::from_json(
                        args.sample.clone(), 
                        path
                    )?,
                    None => GpConfig::new(
                        args.sample.clone(), 
                        args.controls.clone(), 
                        args.tags.clone(),
                        PrevalenceContaminationConfig::gp_default()
                    )
                };
                            
                // Prime the LLM with the background explanation
                agent.prime_llm(
                    &DataBackground::CerebroFilter.get_default()
                ).await?;
                            
                // Run the diagnostic agent, for synthesis of metagenomic pathogen 
                // detection, host aneuploidy signature, and patient clinical notes
                let diagnostic_result = agent.run(
                    &mut gp_config, 
                    &args.clinical_context.text()
                ).await?;

                diagnostic_result.to_json(&args.diagnostic_log)?;

                if let Some(state_log) = &args.state_log {
                    agent.state.to_json(state_log)?;
                }

            }
            

        },
    }

    Ok(())

}
   