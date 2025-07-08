
use std::{fs::create_dir_all, process::exit, sync::Arc};

use meta_gpt::{gpt::{AgentBenchmark, DiagnosticAgent, DiagnosticResult, GpuBenchmark, PrefetchData, SampleContext, TaskConfig},  utils::get_config};

#[cfg(feature = "local")]
use meta_gpt::text::{GeneratorConfig, TextGenerator};

use cerebro_model::api::cerebro::{model::Cerebro, schema::PostFilterConfig};
use cerebro_pipeline::{modules::{pathogen::PathogenDetection, quality::{QualityControl, QualityControlSummary, PositiveControlConfig, PositiveControlSummaryBuilder}}, utils::{get_file_component, FileComponent}};
use clap::Parser;
use cerebro_ciqa::{error::CiqaError, plate::{aggregate_reference_plates, get_diagnostic_stats, load_diagnostic_stats_from_files, plot_plate, plot_qc_summary_matrix, plot_stripplot, DiagnosticData, MissingOrthogonal, Palette, ReferencePlate, SampleReference, SampleType}, stats::mcnemar_from_reviews, terminal::{App, Commands}, utils::{init_logger, write_tsv}};
use cerebro_client::{client::CerebroClient, error::HttpClientError};
use plotters::prelude::SVGBackend;
use plotters_bitmap::BitMapBackend;
use serde::Serialize;
use std::sync::atomic::{AtomicBool, Ordering};

#[tokio::main]
async fn main() -> anyhow::Result<(), anyhow::Error> {
    
    init_logger();

    let cli = App::parse();

    match cli.command {
        Commands::PlotPlate( args ) => {
            plot_plate()?
        },
        Commands::PlotPrefetch( args ) => {
        
            let prefetch = PrefetchData::from_json(args.prefetch)?;
            prefetch.plot_domain_counts_svg(&args.output)?;
        },
        Commands::PlotQc( args ) => {

            let mut qc_summaries = Vec::new();

            if !args.cerebro.is_empty() {
                log::info!(
                    "Reading database models and saving quality control summaries to: {}", 
                    args.outdir.display()
                );

                create_dir_all(&args.outdir)?;

                for file in args.cerebro {
                    let model = Cerebro::from_json(&file)?;
    
                    let qc_summary = model.quality.summary_illumina_pe();
    
                    qc_summary.to_json(&args.outdir.join(
                        format!("qc_{}", get_file_component(
                            &file, 
                            FileComponent::FileName
                        )?)
                    ))?;
    
                    qc_summaries.push(qc_summary);
                }
            } else if !args.summaries.is_empty() {
                log::info!("Reading quality control summaries");

                for file in args.summaries {
                    let qc_summary = QualityControlSummary::from_json(&file)?;
                    qc_summaries.push(qc_summary);
                }
            } else if !args.quality_control.is_empty() {
                log::info!(
                    "Reading quality controle module files and saving quality control summaries to: {}", 
                    args.outdir.display()
                );

                create_dir_all(&args.outdir)?;

                for file in args.quality_control {
                    let model = QualityControl::from_json(&file)?;
    
                    let qc_summary = model.summary_illumina_pe();
    
                    qc_summary.to_json(&args.outdir.join(
                        format!("qc_{}", get_file_component(
                            &file, 
                            FileComponent::FileName
                        )?)
                    ))?;
    
                    qc_summaries.push(qc_summary);
                }

            } else {
                log::error!("Either Cerebro model files (.json, --cerebro, -c), quality control files (.json, --quality-control, -q) or quality control summary files (.json, --summaries, -s) must be provided!");
                exit(1);
            };
    
            plot_qc_summary_matrix(
                &qc_summaries, 
                None, 
                &args.output, 
                cerebro_ciqa::plate::CellShape::Circle, 
                args.width, 
                args.height, 
                args.title.as_deref()
            )?;


            // Positive controls:

            if !args.pathogen_detection.is_empty() {
                log::info!(
                    "Reading pathogen detection models and saving positive control summaries to: {}", 
                    args.outdir.display()
                );

                for file in args.pathogen_detection {
                    let model = PathogenDetection::from_json(&file)?;
                    let summary = PositiveControlSummaryBuilder::new(
                        &model, 
                        &PositiveControlConfig::metagp(),
                        &model.id
                    ).build();
                }

            }

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
            let highlight = Palette::paquin();

            let col1 = palette.colors.get(3);
            let col2 = palette.colors.get(6);
            let col3 = highlight.colors.get(3);

            if ext == "svg" {
                plot_stripplot(
                    SVGBackend::new(&args.output, (args.width, args.height)),
                    &data, 
                    args.mode, 
                    args.ref1, 
                    args.ref2,
                    col1,
                    col2,
                    col3,
                    args.ci,
                    args.stats.clone(),
                    None,
                    args.boxplot,
                    args.barplot,
                    args.y_labels
                )?;
            } else {
                plot_stripplot(
                    BitMapBackend::new(&args.output, (args.width, args.height)),
                    &data, 
                    args.mode, 
                    args.ref1, 
                    args.ref2,
                    col1,
                    col2,
                    col3,
                    args.ci,
                    args.stats.clone(),
                    None,
                    args.boxplot,
                    args.barplot,
                    args.y_labels
                )?;
            };

        },
        Commands::Prefetch( args ) => {

            // Deconvolute the blocking task spawn in this async context, if not (and if this runs in --release, this may be
            // the reason for the ongoing dropouts in requests to the server)

            // In Commands::Prefetch, we are blocking IO and sync operations like create_dir_all, ReferencePlate::new, 
            // and loops with potentially heavy filesystem or network IO inside an async function.
            // Even though prefetch() appears async-safe, ReferencePlate creation, file checking, and even DiagnosticAgent 
            // usage might contain sync-blocking behavior, particularly via filesystem or network access.

            tokio::task::spawn_blocking(move || -> Result<(), anyhow::Error> {

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

                    if let Some(ref subset) = args.samples {
                        if !subset.contains(sample_id) {
                            log::warn!("Sample subset provided - skipping {}", sample_id);
                            continue;
                        }
                    }

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

                        match DiagnosticAgent::new(Some(api_client.clone()), TaskConfig::Default)?
                            .prefetch(&data_file, &gp_config) {
                                Ok(_) => {},
                                Err(err) => match err {
                                    meta_gpt::error::GptError::CerebroClientError(HttpClientError::ResponseFailure(status))
                                     => {
                                        // Sample not found warning but pass by skipping
                                        log::warn!("Skipping sample due to failure to prefetch threshold data: {sample_id}")
                                    }
                            
                                    // re-propagate
                                    other => return Err(other.into()),
                                },
                            }

                    } else {
                        log::info!("File '{}' exists and force not enabled - skipping file", data_file.display());
                    }
                }

                Ok(())
            })
            .await??;
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

        #[cfg(feature = "local")]
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
            
            // How many GPUs (default to 1 if not set)
            let num_gpus = args.num_gpu;
            
            // Collect samples and compute batch size
            let samples: Vec<_> = plate.samples.clone();
            let batch_size = (samples.len() + num_gpus - 1) / num_gpus;
            
            // Share args + plate between threads
            let args = Arc::new(args);
            let plate = Arc::new(plate);
            let nvml = Arc::new(nvml_wrapper::Nvml::init()?);
            
            // Spawn one thread per GPU
            let mut handles = Vec::with_capacity(num_gpus);
            
            for gpu_id in 0..num_gpus {

                let args= Arc::clone(&args);
                let plate  = Arc::clone(&plate);
                let nvml_clone = Arc::clone(&nvml);

                let start  = gpu_id * batch_size;
                let end    = (start + batch_size).min(samples.len());
                let batch_samples = samples[start..end].to_vec();
                
                handles.push(std::thread::spawn(move || -> Result<(), CiqaError> {

                    let keep_polling = Arc::new(AtomicBool::new(true));
                    let polling_flag = Arc::clone(&keep_polling);

                    let bench_handle = std::thread::spawn(move || -> Result<u64, CiqaError> {

                        let nvml_device = nvml_clone.device_by_index(gpu_id as u32)?;
                        let mut peak = 0;

                        while polling_flag.load(Ordering::Relaxed) {
                            let usage = nvml_device.memory_info()?.used;
                            peak = peak.max(usage);
                            std::thread::sleep(std::time::Duration::from_millis(100));
                        }

                        Ok(peak)
                    });

                    // state logs directory
                    let state_dir = args.outdir.join("state_logs");
                    let gpu_bench_file = state_dir.join(format!("gpu{gpu_id}.bench.json"));

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

                    for sample_id in batch_samples {
                        
                        // rebuild paths

                        use meta_gpt::gpt::ClinicalContext;
                        let data_file   = args.prefetch.join(format!("{sample_id}.prefetch.json"));

                        let result_file = args.outdir.join(format!("{sample_id}.model.json"));
                        let state_file  = state_dir.join(format!("{sample_id}.state.json"));
                        let bench_file  = state_dir.join(format!("{sample_id}.bench.json"));
                        
                        if !data_file.exists() {
                            log::info!("no prefetch for {}, skipping", sample_id);
                            continue;
                        }
                        
                        // load prefetch
                        let prefetch_data = PrefetchData::from_json(&data_file)?;
                        log::info!("prefetch.config = {:#?}", prefetch_data.config);

                        // instantiate agent
                        let mut agent = DiagnosticAgent::new(None, args.task_config.clone())?;
                        
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
                        
                        let clinical_notes = sample_ref.and_then(|s: SampleReference| {
                            if let Some(clinical) = s.clinical {
                                Some(ClinicalContext::Custom(clinical))
                            } else {
                                None
                            }
                        });
                        
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

                        let start = std::time::Instant::now();

                        // run agent
                        let result = agent.run_local(
                            &mut text_generator,
                            if args.sample_context.unwrap_or(false) { Some(sample_context) } else { None },
                            if args.clinical_notes { clinical_notes } else { None },
                            args.assay_context.clone(),
                            args.agent_primer.clone(),
                            &prefetch_data.config.clone(),
                            Some(prefetch_data),
                            post_filter,
                            args.disable_thinking
                        )?;

                        let elapsed = start.elapsed().as_secs_f32();

                        let bench = AgentBenchmark {
                            seconds: elapsed
                        };

                        // write out
                        result.to_json(&result_file)?;
                        agent.state.to_json(&state_file)?;
                        bench.to_json(&bench_file)?;

                        log::info!("finished {}", sample_id);
                    }

                    // Stop polling
                    keep_polling.store(false, Ordering::Relaxed);

                    let peak_vram = match bench_handle.join() {
                        // Thread ran to completion, giving a `Result<(), GptError>`
                        Ok(inner) => {
                            let peak_vram = match inner {
                                Ok(peak_vram) => peak_vram,
                                Err(gpt_err) => {
                                    log::error!("GPT error in GPU memory thread: {}", gpt_err.to_string());
                                    0
                                }
                            };
                            peak_vram
                        }
                        // Thread actually panicked:
                        Err(panic_payload) => {
                            log::error!("GPU memory thread panicked: {:?}", panic_payload);
                            0
                        }
                    };

                    GpuBenchmark {  peak_vram }.to_json(&gpu_bench_file)?;

                    Ok(())
                }));
            }
            
            // Wait for all GPUs to finish
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

        Commands::Mcnemar( args ) => {

            let data_a = DiagnosticData::from_json(args.review_a)?;
            let data_b = DiagnosticData::from_json(args.review_b)?;

            let consensus_reviews_a = data_a.get_consensus_reviews_filtered();
            let consensus_reviews_b = data_b.get_consensus_reviews_filtered();

            match (consensus_reviews_a, consensus_reviews_b) {
                (Some(reviews_a), Some(reviews_b)) => {
                    let test_result = mcnemar_from_reviews(&reviews_a, &reviews_b);
                    log::info!("{test_result:#?}")
                },
                _ => {
                    log::error!("Could not find consensus reviews in provided diagnostic review data!")
                }
            }

        }
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
