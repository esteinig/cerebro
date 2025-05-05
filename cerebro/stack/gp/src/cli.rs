
use std::path::PathBuf;

use cerebro_model::api::cerebro::schema::{MetaGpConfig, PrevalenceOutliers};
use cerebro_pipeline::taxa::filter::PrevalenceContaminationConfig;
use clap::Parser;
use cerebro_gp::{error::GptError, gpt::{draw_consensus_tree, AssayContext, DiagnosticAgent, PrefetchData}, terminal::{App, Commands}, utils::init_logger};
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
                args.contam_history,
                None
            )?;

            log::info!("{:#?}", gp_config);

            DiagnosticAgent::new(Some(api_client))
                .await?
                .prefetch(&args.output, &gp_config)
                .await?
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

                        
            let mut agent = DiagnosticAgent::new(api_client).await?;


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
                args.contam_history,
                args.prefetch.clone()
            )?;

            let result = agent.run_local(
                &mut generator,
                args.sample_context.clone(),
                args.clinical_notes.clone(),
                args.assay_context.clone(),
                &gp_config, 
                prefetch
            ).await?;

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


fn get_config(
    json: &Option<PathBuf>,
    sample: Option<String>,
    controls: &Option<Vec<String>>,
    tags: &Option<Vec<String>>,
    ignore_taxstr: &Option<Vec<String>>,
    prevalence_outliers: Option<bool>,
    prefetch: Option<PathBuf>
) -> Result<(MetaGpConfig, Option<PrefetchData>), GptError> {

    // Parse generative practitioner configuration
    match &json {
        Some(path) => {
            // Config from JSON file:
            Ok((
                MetaGpConfig::with_json(
                    sample.ok_or(GptError::SampleIdentifierMissing)?, 
                    path,
                    ignore_taxstr.clone(),
                    match prevalence_outliers { 
                        Some(true) => Some(PrevalenceOutliers::default()),
                        Some(false) => Some(PrevalenceOutliers::disabled()),
                        None => None
                    }
                )?,
                None
            ))
        },
        None => {
            match &prefetch {
                Some(path) => {
                    let data = PrefetchData::from_json(path)?;
                    Ok((data.config.clone(), Some(data)))
                },
                None => {
                    Ok((
                        MetaGpConfig::new(
                            sample.ok_or(GptError::SampleIdentifierMissing)?, 
                            controls.clone(), 
                            tags.clone(),
                            ignore_taxstr.clone(),
                            PrevalenceContaminationConfig::gp_default(),
                            match prevalence_outliers { 
                                Some(true) => Some(PrevalenceOutliers::default()),
                                Some(false) => Some(PrevalenceOutliers::disabled()),
                                None => None
                            }
                        ),
                        None
                    ))
                }
            }
        }
    }
}