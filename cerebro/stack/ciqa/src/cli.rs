
use std::fs::create_dir_all;

use cerebro_gp::gpt::{ClinicalContext, DataBackground, DiagnosticAgent};
use cerebro_model::api::cerebro::schema::{CerebroIdentifierSchema, GpConfig};
use cerebro_pipeline::taxa::filter::PrevalenceContaminationConfig;
use clap::Parser;
use cerebro_ciqa::{config::EvaluationConfig, plate::{plot_plate, DiagnosticOutcome, FromSampleType, MissingOrthogonal, ReferencePlate}, terminal::{App, Commands}, utils::init_logger};
use cerebro_client::client::CerebroClient;

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

            api_client.get_taxa(
                &request_schema, 
                &eval_config.filter, 
                &mut eval_config.prevalence,
                args.contam_history
            )?;
        },
        Commands::Plate( args ) => {

            plot_plate()?
        },
        Commands::Review( args ) => {

            let mut plate_reference = ReferencePlate::new(
                &args.plate, 
                args.review.as_deref(),
                args.missing_orthogonal.clone(),
                args.diagnostic_agent
            )?;

            if let Some(sample_id) = &args.set_none {
                plate_reference.set_none(sample_id)?;
            }

            if let Some(sample_type) = &args.sample_type {
                plate_reference.subset_sample_type(sample_type.clone())?;
            }

            let diagnostic_review = plate_reference.compute_diagnostic_review()?;

            for dr in &diagnostic_review {
                log::info!("{} => {}", dr.sample_id, dr.outcome.colored())
            }

            let stats = plate_reference.compute_statistics(diagnostic_review);

            log::info!("{:#?}", stats);
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


            for sample_id in &plate.samples {
                
                let json_file = &args.outdir.join(
                    format!("{sample_id}.{}.json", args.model)
                );

                if args.force || !json_file.exists() {

                    let mut agent = DiagnosticAgent::new(
                        api_client.clone(), 
                        args.model.clone(),
                        args.diagnostic_memory,
                        args.contam_history
                    ).await?;
        
                    
                    // Prime the LLM with the background explanation
                    agent.prime_llm(
                        &DataBackground::CerebroFilter.get_default()
                    ).await?;
    
                    let sample_reference = match plate.get_sample_reference(&sample_id) {
                        Some(sample_reference) => sample_reference,
                        None => { 
                            log::warn!("Failed to find plate reference for sample: {sample_id}");
                            continue;
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
                                Some(note.replace("ignore_taxstr:", "").split(",").map(|s| s.trim().to_string()).collect::<Vec<_>>())
                            } else {
                                None
                            }
                        },
                        None => None
                    };
                    
                    let mut gp_config = GpConfig::from_reference_plate(
                        sample_id, 
                        &plate.negative_controls, 
                        &tags,
                        ignore_taxstr,
                        PrevalenceContaminationConfig::gp_default()
                    );
    
                    let diagnostic_result = agent.run(
                        &mut gp_config, 
                        &ClinicalContext::from_sample_type(sample_reference.sample_type).text()
                    ).await?;
    
                    diagnostic_result.to_json(&json_file)?;
                } else {
                    log::info!("File '{}' exists and force not enabled - skipping file", json_file.display())
                }

                            
            }
                        
        } 
    }

    Ok(())

}
   