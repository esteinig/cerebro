
use std::path::PathBuf;

use cerebro_model::api::cerebro::schema::{MetaGpConfig, PostFilterConfig, PrevalenceOutliers};
use cerebro_pipeline::taxa::filter::PrevalenceContaminationConfig;
use clap::Parser;
use cerebro_gp::{error::GptError, gpt::{draw_consensus_tree, AssayContext, DiagnosticAgent, PrefetchData}, terminal::{App, Commands}, utils::{get_config, init_logger}};
use cerebro_client::client::CerebroClient;

#[cfg(feature = "local")]
use cerebro_gp::text::{TextGenerator, GeneratorConfig};

#[tokio::main]
async fn main() -> anyhow::Result<(), anyhow::Error> {
    
    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::DiagnoseApi( args ) => {
            log::warn!("Not implemented yet")
        }
        Commands::PrefetchTiered( args ) => {
           
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

            let (gp_config, _) = get_config(
                &args.json, 
                Some(args.sample.clone()), 
                &args.controls, 
                &args.tags, 
                &args.ignore_taxstr,
                args.prevalence_outliers,
                None,   // not relevant
            )?;

            log::info!("{:#?}", gp_config);

            DiagnosticAgent::new(Some(api_client))?
                .prefetch(&args.output, &gp_config)?
        }
        #[cfg(feature = "local")]
        Commands::DiagnoseLocal( args ) => {

            let api_client = match args.prefetch {
                Some(_) => None,
                None => {
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

                    Some(api_client)
                },
            };

                        
            let mut agent = DiagnosticAgent::new(api_client)?;


            let mut generator = TextGenerator::new(
                GeneratorConfig::with_default(
                    args.model.clone(),
                    args.model_dir.clone(),
                    args.sample_len,
                    args.temperature,
                    args.gpu
                )
            )?;

            let (gp_config, prefetch) = get_config(
                &args.json, 
                args.sample.clone(), 
                &args.controls, 
                &args.tags, 
                &args.ignore_taxstr,
                args.prevalence_outliers,
                args.prefetch.clone()
            )?;

            
            let post_filter = if args.post_filter.unwrap_or(false) { 
                let collapse_variants = match args.collapse_variants {
                    Some(value) => value,
                    None => false,
                };
                let species_domains = match &args.species_domains {
                    Some(domains) => domains.to_vec(),
                    None => vec!["Archaea".to_string(), "Bacteria".to_string(), "Eukaryota".to_string()]
                };
                let exclude_phage = match args.exclude_phage {
                    Some(value) => value,
                    None => false,
                };
                Some(PostFilterConfig::with_default(
                    collapse_variants,
                    args.min_species > 0,
                    args.min_species,
                    species_domains,
                    exclude_phage
                )) 
            } else { 
                None 
            };

            let result = agent.run_local(
                &mut generator,
                args.sample_context.clone(),
                args.clinical_notes.clone(),
                args.assay_context.clone(),
                &gp_config, 
                prefetch,
                post_filter
            )?;

            result.to_json(&args.diagnostic_log)?;

            log::info!("{:#?}", result);

            if let Some(state_log) = &args.state_log {
                agent.state.to_json(state_log)?;
            }
            

        },
        #[cfg(feature = "local")]
        Commands::Generate( args ) => {
            
            let mut generator = TextGenerator::new(
                GeneratorConfig::from_args(&args)
            )?;

            generator.run(&args.prompt)?;

        }
    }

    Ok(())

}

