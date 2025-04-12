
use std::fs::create_dir_all;

use cerebro_gp::gpt::{ClinicalContext, DataBackground, DiagnosticAgent};
use cerebro_model::api::cerebro::schema::{CerebroIdentifierSchema, MetaGpConfig};
use cerebro_pipeline::taxa::filter::PrevalenceContaminationConfig;
use clap::Parser;
use cerebro_ciqa::{config::EvaluationConfig, plate::{aggregate_reference_plates, get_diagnostic_stats, load_diagnostic_stats_from_files, plot_plate, plot_stripplot, read_all_sample_reviews, DiagnosticStatsVecExt, FromSampleType, MissingOrthogonal, Palette, ReferencePlate}, terminal::{App, Commands}, utils::{get_file_component, init_logger, FileComponent}};
use cerebro_client::client::CerebroClient;
use plotters::prelude::SVGBackend;
use plotters_bitmap::BitMapBackend;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use tokio::runtime::Runtime;

#[tokio::main]
async fn main() -> anyhow::Result<(), anyhow::Error> {
    
    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Evaluate( args ) => {

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
            
            let mut eval_config = match &args.json {
                Some(path) => EvaluationConfig::from_json(path)?,
                None => EvaluationConfig::from_evaluate_args(args)
            };
            
            let request_schema = CerebroIdentifierSchema::new(
                args.sample.clone(),
                (!eval_config.controls.is_empty()).then(|| eval_config.controls.clone()),
                (!eval_config.tags.is_empty()).then(|| eval_config.tags.clone()),
            );

            let (taxa, contam_taxa) = api_client.get_taxa(
                &request_schema, 
                &eval_config.filter, 
                &mut eval_config.prevalence,
                args.contam_history
            )?;
            
            for taxon in contam_taxa {
                log::info!("Prevalence contamination taxon: {:#?}", taxon);
            }

            for taxon in taxa {
                log::info!("Detected taxon: {:#?}", taxon);
            }

        },
        Commands::PlotPlate( args ) => {
            plot_plate()?
        },
        Commands::PlotReview( args ) => {
            
            let data = load_diagnostic_stats_from_files(
                args.stats.clone()
            )?;


            let ext = args.output
                .extension()
                .unwrap_or_default()
                .to_string_lossy()
                .to_lowercase();

            let palette = Palette::hiroshige();

            if ext == "svg" {
                plot_stripplot(
                    SVGBackend::new(&args.output, (args.width, args.height)),
                    &data, 
                    args.mode, 
                    args.ref1, 
                    args.ref2,
                    palette.colors.get(4),
                    palette.colors.get(5)
                )?;
            } else {
                plot_stripplot(
                    BitMapBackend::new(&args.output, (args.width, args.height)),
                    &data, 
                    args.mode, 
                    args.ref1, 
                    args.ref2,
                    palette.colors.get(4),
                    palette.colors.get(5)
                )?;
            };

        },
        Commands::Review( args ) => {

            let mut review_data = Vec::new();
            let mut reference_plates = Vec::new();

            for review_path in &args.review {
                
                let review_name = get_file_component(
                    review_path, 
                    FileComponent::FileStem
                )?;
                log::info!("Review: {}", review_name);

                let mut reference_plate = ReferencePlate::new(
                    &args.plate, 
                    Some(&review_path),
                    args.missing_orthogonal.clone(),
                    args.diagnostic_agent
                )?;

                let stats = get_diagnostic_stats(
                    args,
                    &mut reference_plate, 
                    &review_name
                )?;

                let stats_percent = stats.percent();
                log::info!("{stats_percent:#?}");

                review_data.push(stats_percent);
                reference_plates.push(reference_plate);
            }
            
            log::info!("Consensus (majority vote) review tagged 'consensus'");

            let mut consensus_plate = aggregate_reference_plates(reference_plates);

            let stats = get_diagnostic_stats(
                args,
                &mut consensus_plate, 
                &"consensus"
            )?;

            let stats_percent = stats.percent();
            log::info!("{stats_percent:#?}");

            review_data.push(stats_percent);
            review_data.to_json(&args.output)?;

        },
        Commands::Diagnose( args ) => {

            create_dir_all(&args.outdir)?;

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


            let plate = ReferencePlate::new(
                &args.plate, 
                None,
                MissingOrthogonal::Indeterminate,
                false
            )?;

            let pool = ThreadPoolBuilder::new()
                .num_threads(args.threads)
                .build()
                .expect("Failed to create the Rayon thread pool");

            pool.install(|| {
                plate.samples.par_iter().for_each(|sample_id| {
                    // Create a new runtime in each thread
                    let rt = Runtime::new().expect("Failed to create Tokio runtime");

                    // Run the async operations inside the synchronous Rayon thread
                    let result: Result<(), anyhow::Error> = rt.block_on(async {

                        let json_file = args.outdir.join(format!("{sample_id}.{}.json", args.model));

                        if args.force || !json_file.exists() {
                            let mut agent = DiagnosticAgent::new(
                                api_client.clone(), 
                                args.model.clone(),
                                args.diagnostic_memory,
                                args.contam_history
                            ).await?;

                            // Prime the LLM with the background explanation
                            agent.prime_llm(&DataBackground::CerebroFilter.get_default()).await?;

                            let sample_reference = match plate.get_sample_reference(&sample_id) {
                                Some(sample_reference) => sample_reference,
                                None => { 
                                    log::warn!("Failed to find plate reference for sample: {sample_id}");
                                    return Ok(());
                                }
                            };

                            log::info!("Starting generative practitioner diagnostic for '{sample_id}'");

                            let tags = match sample_reference.note {
                                Some(ref note) => {
                                    if note.contains("ignore_rna") {
                                        vec![String::from("DNA")]
                                    } else if note.contains("ignore_dna") {
                                        vec![String::from("RNA")]
                                    } else {
                                        vec![String::from("DNA"), String::from("RNA")]
                                    }
                                },
                                None => vec![String::from("DNA"), String::from("RNA")]
                            };

                            let ignore_taxstr = match sample_reference.note {
                                Some(ref note) => {
                                    if note.contains("ignore_taxstr") {
                                        Some(note.replace("ignore_taxstr:", "")
                                            .split(',')
                                            .map(|s| s.trim().to_string())
                                            .collect::<Vec<_>>())
                                    } else {
                                        None
                                    }
                                },
                                None => None
                            };

                            let mut meta_gp_config = MetaGpConfig::from_reference_plate(
                                sample_id, 
                                &plate.negative_controls, 
                                &tags,
                                ignore_taxstr,
                                PrevalenceContaminationConfig::gp_default()
                            );

                            let diagnostic_result = agent.run(
                                &mut meta_gp_config, 
                                &ClinicalContext::from_sample_type(sample_reference.sample_type).text()
                            ).await?;

                            diagnostic_result.to_json(&json_file)?;
                        } else {
                            log::info!("File '{}' exists and force not enabled - skipping file", json_file.display());
                        }
                        Ok(())
                    });
                    if let Err(e) = result {
                        log::error!("Error processing sample {}: {:?}", sample_id, e);
                    }
                });
            });
                        
        } 
    }

    Ok(())

}
   