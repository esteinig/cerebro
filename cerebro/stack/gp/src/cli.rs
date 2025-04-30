
use cerebro_model::api::cerebro::schema::MetaGpConfig;
use cerebro_pipeline::taxa::filter::PrevalenceContaminationConfig;
use clap::Parser;
use cerebro_gp::{gpt::{draw_consensus_tree, DataBackground, DiagnosticAgent}, terminal::{App, Commands}, utils::init_logger};
use cerebro_client::client::CerebroClient;

#[cfg(feature = "local")]
use cerebro_gp::llama::run_llama;
#[cfg(feature = "local")]
use cerebro_gp::quantized::run_quantized;
#[cfg(feature = "local")]
use cerebro_gp::qwen::run_qwen;
#[cfg(feature = "local")]
use cerebro_gp::text::TextGenerator;

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
                args.model.clone(),
                args.diagnostic_memory,
                args.contam_history
            ).await?;

            // Optionally, print the knowledge graph
            agent.print_decision_tree(&agent.tree, "start");

            draw_consensus_tree(&DiagnosticAgent::graph(&agent.tree)?, &vec![], "tree.svg", 1000, 1000)?;

            if !args.dry_run {

                // Parse generative practicioner configuration
                let mut gp_config = match &args.json {
                    Some(path) => {
                        MetaGpConfig::from_json(
                            args.sample.clone(), 
                            path,
                            args.ignore_taxstr.clone()
                        )?
                    },
                    None => MetaGpConfig::new(
                        args.sample.clone(), 
                        args.controls.clone(), 
                        args.tags.clone(),
                        args.ignore_taxstr.clone(),
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
        #[cfg(feature = "local")]
        Commands::Generate( args ) => {
            
            let generator = TextGenerator::new(
                args.model, 
                &args.model_dir, 
                args.force_download
            );

            generator.run(
                &args.prompt,
                args.raw_prompt,
                args.sample_len,
                args.split_prompt,
                args.clean,
                args.temperature,
                args.seed,
                args.top_k,
                args.top_p,
                args.repeat_penalty,
                args.repeat_last_n,
                args.log_info
            )?;

        }
        #[cfg(feature = "local")]
        Commands::Llama( args ) => {
            run_llama(args.clone())?
        },
        #[cfg(feature = "local")]
        Commands::Quantized( args ) => {
            run_quantized(args.clone())?
        },
        #[cfg(feature = "local")]
        Commands::Qwen( args ) => {
            run_qwen(args.clone())?
        }
    }

    Ok(())

}
   