
use std::{fs::create_dir_all, path::PathBuf, sync::Arc};

use cerebro_gp::{error::GptError, gpt::{AssayContext, DiagnosticAgent, DiagnosticResult, PrefetchData, SampleContext}, text::{GeneratorConfig, TextGenerator}, utils::get_config};
use cerebro_model::api::cerebro::schema::{CerebroIdentifierSchema, MetaGpConfig, PostFilterConfig, PrevalenceOutliers};
use cerebro_pipeline::taxa::filter::PrevalenceContaminationConfig;
use clap::Parser;
use cerebro_ciqa::{config::EvaluationConfig, error::CiqaError, plate::{aggregate_reference_plates, get_diagnostic_stats, load_diagnostic_stats_from_files, plot_diagnostic_matrix, plot_plate, plot_stripplot, CellShape, DiagnosticData, DiagnosticReview, DiagnosticStatsVecExt, DiagnosticSummary, FromSampleType, MissingOrthogonal, Palette, ReferencePlate, SampleReference, SampleType}, terminal::{App, Commands}, utils::{get_file_component, init_logger, write_tsv, FileComponent}};
use cerebro_client::client::CerebroClient;
use plotters::prelude::SVGBackend;
use plotters_bitmap::BitMapBackend;
use serde::Serialize;

#[tokio::main]
async fn main() -> anyhow::Result<(), anyhow::Error> {
    
    init_logger();

    let cli = App::parse();

    match cli.command {
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
            let highlight = Palette::cassatt2();

            if ext == "svg" {
                plot_stripplot(
                    SVGBackend::new(&args.output, (args.width, args.height)),
                    &data, 
                    args.mode, 
                    args.ref1, 
                    args.ref2,
                    palette.colors.get(3),
                    palette.colors.get(6),
                    highlight.colors.get(3),
                    args.ci
                )?;
            } else {
                plot_stripplot(
                    BitMapBackend::new(&args.output, (args.width, args.height)),
                    &data, 
                    args.mode, 
                    args.ref1, 
                    args.ref2,
                    palette.colors.get(3),
                    palette.colors.get(6),
                    highlight.colors.get(3),
                    args.ci
                )?;
            };

        },
        Commands::Prefetch( args ) => {


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

            create_dir_all(&args.outdir)?;

            let plate = ReferencePlate::new(
                &args.plate, 
                None,
                MissingOrthogonal::Indeterminate,
                false
            )?;

            for sample_id in &plate.samples {

                let data_file = args.outdir.join(format!("{sample_id}.prefetch.json"));

                if args.force || !data_file.exists() {

                    let sample_reference = match plate.get_sample_reference(&sample_id) {
                        Some(sample_reference) => sample_reference,
                        None => { 
                            log::warn!("Failed to find plate reference for sample: {sample_id}");
                            break
                        }
                    };

                    let (tags, ignore_taxstr) = get_note_instructions(&sample_reference);

                    let (gp_config, _) = get_config(
                        &None,
                        Some(sample_reference.sample_id), 
                        &Some(plate.negative_controls.clone()), 
                        &Some(tags), 
                        &ignore_taxstr,
                        args.prevalence_outliers,
                        None
                    )?;

                    log::info!("{:#?}", gp_config);

                    DiagnosticAgent::new(Some(api_client.clone()))?
                        .prefetch(&data_file, &gp_config)?;

                } else {
                    log::info!("File '{}' exists and force not enabled - skipping file", data_file.display());
                }
            }
        }
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
                    &args,
                    &mut reference_plate, 
                    &review_name
                )?;

                let stats_percent = stats.percent();
                log::info!("{stats_percent}");

                review_data.push(stats_percent);
                reference_plates.push(reference_plate);
            }

            log::info!("Consensus (majority vote) review tagged 'consensus'");
            
            let mut consensus_plate = aggregate_reference_plates(reference_plates);

            let consensus_stats = get_diagnostic_stats(
                &args,
                &mut consensus_plate, 
                &"consensus"
            )?;

            let consensus_stats_percent = consensus_stats.percent();
            log::info!("{consensus_stats_percent}");

            review_data.push(consensus_stats_percent);

            let data = DiagnosticData::from(review_data);

            data.to_json(&args.output)?;

            if let Some(plot_path) = &args.plot {
                data.plot_summary(
                    plot_path, 
                    args.title.as_deref(), 
                    args.width, 
                    args.height, 
                    args.reference.clone()
                )?;
            }
        },
        Commands::DiagnoseLocal( args ) => {

            std::fs::create_dir_all(&args.outdir)?;

            let state_dir = args.outdir.join("state_logs");
            std::fs::create_dir_all(&state_dir)?;

            let plate = ReferencePlate::new(
                &args.plate,
                None,
                MissingOrthogonal::Indeterminate,
                false,
            )?;
            
            // 2) How many GPUs (default to 1 if not set)
            let num_gpus = args.num_gpu;
            
            // 3) Collect samples and compute batch size
            let samples: Vec<_> = plate.samples.clone();
            let batch_size = (samples.len() + num_gpus - 1) / num_gpus;
            
            // 4) Share args + plate between threads
            let args = Arc::new(args);
            let plate = Arc::new(plate);
            
            // 5) Spawn one thread per GPU
            let mut handles = Vec::with_capacity(num_gpus);
            
            for gpu_id in 0..num_gpus {

                let args   = Arc::clone(&args);
                let plate  = Arc::clone(&plate);

                let start  = gpu_id * batch_size;
                let end    = (start + batch_size).min(samples.len());
                let batch_samples = samples[start..end].to_vec();
                
                handles.push(std::thread::spawn(move || -> Result<(), GptError> {

                    // state logs directory
                    let state_dir = args.outdir.join("state_logs");

                    // load inference model and weights for this batch on GPU
                    let mut text_generator = TextGenerator::new(
                        GeneratorConfig::with_default(
                            args.model.clone(),
                            args.model_dir.clone(),
                            args.sample_len,
                            args.temperature,
                            gpu_id,
                        )
                    )?;

                    // post filter config log


                    for sample_id in batch_samples {
                        
                        // rebuild paths
                        let data_file   = args.prefetch.join(format!("{sample_id}.prefetch.json"));

                        let result_file = args.outdir.join(format!("{sample_id}.model.json"));
                        let state_file  = state_dir.join(format!("{sample_id}.state.json"));
                        
                        if !data_file.exists() {
                            log::info!("no prefetch for {}, skipping", sample_id);
                            continue;
                        }
                        
                        // load prefetch
                        let prefetch_data = PrefetchData::from_json(&data_file)?;
                        log::info!("prefetch.config = {:#?}", prefetch_data.config);

                        // instantiate agent
                        let mut agent = DiagnosticAgent::new(None)?;
                        
                        // sample context logic
                        let sample_ref = plate.get_sample_reference(&sample_id);

                        let sample_context = sample_ref
                            .as_ref()
                            .map(|s| match s.sample_type {
                                SampleType::Csf => SampleContext::Csf,
                                SampleType::Eye => SampleContext::Eye,
                                _               => SampleContext::None,
                            })
                            .unwrap_or(SampleContext::None);
                        
                        let clinical_notes = sample_ref.and_then(|s| s.clinical.clone());
                        
                        let post_filter = if args.post_filter.unwrap_or(false) { 
                            let collapse_variants = match args.collapse_variants {
                                Some(value) => value,
                                None => false,
                            };
                            Some(PostFilterConfig::with_default(
                                collapse_variants,
                                args.min_species > 0,
                                args.min_species
                            )) 
                        } else { 
                            None 
                        };

                        // run agent
                        let result = agent.run_local(
                            &mut text_generator,
                            if args.sample_context.unwrap_or(false) { sample_context } else { SampleContext::None },
                            if args.clinical_notes { clinical_notes } else { None },
                            args.assay_context.clone(),
                            &prefetch_data.config.clone(),
                            Some(prefetch_data),
                            post_filter
                        )?;
                        
                        // write out
                        result.to_json(&result_file)?;
                        agent.state.to_json(&state_file)?;
                        log::info!("finished {}", sample_id);
                    }
                    Ok(())
                }));
            }
            
            // 6) Wait for all GPUs to finish
            for handle in handles {
                match handle.join() {
                    // Thread ran to completion, giving you a `Result<(), GptError>`
                    Ok(inner) => {
                        match inner {
                            Ok(()) => {}
                            Err(gpt_err) => {
                                log::error!("GPT error in thread: {}", gpt_err.to_string());
                            }
                        }
                    }
                    // Thread actually panicked:
                    Err(panic_payload) => {
                        log::error!("Worker thread panicked: {:?}", panic_payload);
                    }
                }
            }                  
        },

        Commands::DebugPathogen( args ) => {

            #[derive(Serialize)]
            struct PathogenRecord {
                filename: String,
                pathogen: String
            }

            let mut records = Vec::new();
            for path in &args.gpt {
                
                let filename = get_file_component(
                    path, 
                    FileComponent::FileStem
                )?;
                log::info!("Diagnostic result: {}", filename);

                let result = DiagnosticResult::from_json(path)?;
                
                records.push(PathogenRecord {
                    filename, 
                    pathogen: result.pathogen.unwrap_or("".to_string())
                });

            }            
            write_tsv(&records, &args.output, false)?;

        }
    }

    Ok(())

}

fn get_note_instructions(sample_reference: &SampleReference) -> (Vec<String>, Option<Vec<String>>) {

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

    (tags, ignore_taxstr)
}