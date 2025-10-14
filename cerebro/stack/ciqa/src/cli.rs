
use std::{collections::{HashMap, HashSet}, fs::create_dir_all, process::exit, sync::Arc};
use cerebro_model::api::cerebro::schema::SampleType;
use cerebro_fs::client::{FileSystemClient, UploadConfig};
use meta_gpt::{gpt::{AgentBenchmark, DiagnosticAgent, DiagnosticResult, GpuBenchmark, SampleContext, TaskConfig, ClinicalContext}};

#[cfg(feature = "local")]
use meta_gpt::text::{GeneratorConfig, TextGenerator};

use cerebro_model::api::{cerebro::{model::Cerebro, schema::{MetaGpConfig, PostFilterConfig, PrefetchData, PrevalenceContaminationConfig, TieredFilterConfig,}}, files::model::FileType};
use cerebro_pipeline::{modules::{pathogen::{PathogenDetection, PathogenDetectionTableRecord}, quality::{write_positive_control_summaries, PositiveControlConfig, PositiveControlSummary, PositiveControlSummaryBuilder, QualityControl, QualityControlSummary}}, utils::{get_file_component, FileComponent}};
use clap::Parser;
use cerebro_ciqa::{error::CiqaError, plate::{aggregate_reference_plates, get_diagnostic_stats, load_diagnostic_stats_from_files, plot_plate, plot_qc_summary_matrix, plot_stripplot, DiagnosticData, MissingOrthogonal, Palette, ReferencePlate, SampleReference}, prefetch::{counts_by_category, counts_by_category_contam, is_missed_detection, positive_candidate_match, reference_names_from_config, MissedDetectionRow, OverallSummary, PerSampleSummary, PrefetchStatus}, stats::{mcnemar_batch_adjust, mcnemar_from_reviews}, terminal::{App, Commands}, utils::{init_logger, write_tsv}};
use cerebro_client::client::CerebroClient;
use plotters::prelude::SVGBackend;
use plotters_bitmap::BitMapBackend;
use serde::Serialize;
use std::sync::atomic::{AtomicBool, Ordering};
use rayon::prelude::*;
use cerebro_ciqa::tui::start_tui;

use std::fs::File;
use std::io::BufReader;

fn main() -> anyhow::Result<(), anyhow::Error> {
    
    init_logger();

    let cli = App::parse();

    match cli.command {
        Commands::PlotPlate( args ) => {
            plot_plate()?
        },
        Commands::PlotPrefetch( args ) => {
            let prefetch = PrefetchData::from_json(args.prefetch)?;
        },
        Commands::WritePlateTable( args ) => {
            let plate = ReferencePlate::from_path(&args.plate)?;
            plate.write_tsv(&args.output, args.species_rank)?;
        },
        Commands::PlotQc( args ) => {


            std::env::set_var("RAYON_NUM_THREADS", args.threads.to_string());

            let summaries: Vec<Result<(QualityControlSummary, PositiveControlSummary), anyhow::Error>> = args.cerebro.par_iter().map(|path| -> anyhow::Result<(QualityControlSummary, PositiveControlSummary), anyhow::Error> {
                log::info!("Reading model file: {}", path.display());

                let model = Cerebro::from_json(path)?;

                let qc_summary = model.quality.summary_illumina_pe();
                let pos_summary = model.summary_positive_control();

                if let Some(outdir) = &args.outdir {
                    create_dir_all(&outdir)?;
                    qc_summary.to_json(&outdir.join(
                        format!("{}", get_file_component(
                            &path, 
                            FileComponent::FileName
                        )?)
                    ))?;
                } 

                Ok((qc_summary, pos_summary))

            }).collect();

            let (qc_summaries, pos_summaries): (Vec<_>, Vec<_>) = summaries.into_iter()
                .filter_map(Result::ok)
                .unzip();

            plot_qc_summary_matrix(
                &qc_summaries, 
                None, 
                &args.output, 
                cerebro_ciqa::plate::CellShape::Circle, 
                args.column_header,
                args.title.as_deref()
            )?;
            
            if let Some(path) = args.positive_controls {
                write_positive_control_summaries(
                    &pos_summaries, 
                    &PositiveControlConfig::metagp(), 
                    path
                )?;
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
        // Deconvolute the blocking task spawn in this async context, if not (and if this runs in --release, this may be
        // the reason for the ongoing dropouts in requests to the server)

        // In Commands::Prefetch, we are blocking IO and sync operations like create_dir_all, ReferencePlate::new, 
        // and loops with potentially heavy filesystem or network IO inside an async function.
        // Even though prefetch() appears async-safe, ReferencePlate creation, file checking, and even DiagnosticAgent 
        // usage might contain sync-blocking behavior, particularly via filesystem or network access.
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

            let plate = ReferencePlate::from_review(
                &args.plate, 
                None,
                MissingOrthogonal::Indeterminate,
                false
            )?;


            let contam_config: PrevalenceContaminationConfig = match args.contamination {
                Some(path) => PrevalenceContaminationConfig::from_json(&path)?,
                None => PrevalenceContaminationConfig::default()
            };

            let prevalence_contamination = if args.disable_filter {
                HashMap::new()
            } else {
                plate.prevalence_contamination(
                    &api_client, 
                    &contam_config
                )?
            };

            let tiered_filter_config = match args.tiered_filter {
                Some(path) => TieredFilterConfig::from_json(&path)?,
                None => {
                    if args.disable_filter {
                        TieredFilterConfig::none()
                    } else {
                        TieredFilterConfig::default(None)
                    }
                }
            };

            std::env::set_var("RAYON_NUM_THREADS", args.threads.to_string());

            let results: Vec<(PerSampleSummary, Option<MissedDetectionRow>)> = plate
                .samples
                .par_iter()
                .filter_map(|sid| {
                    let client = api_client.clone();
                    let outdir = args.outdir.clone();
                    let force = args.force;
                    let subset = args.samples.clone();

                    let contam_config = contam_config.clone();
                    let tiered_filter_config = tiered_filter_config.clone();

                    let sample_id = sid.clone();
                    let negative_controls = plate.negative_controls.clone();
                    let prevalence_contamination = prevalence_contamination.clone();

                    let result: anyhow::Result<Option<(PerSampleSummary, Option<MissedDetectionRow>)>> = (|| {
                        if let Some(ref subset) = subset {
                            if !subset.contains(&sample_id) {
                                log::warn!("Skipping {} (not in subset)", sample_id);
                                return Ok(None);
                            }
                        }

                        let data_file = outdir.join(format!("{sample_id}.prefetch.json"));

                        // Get PrefetchData either by computing or reading
                        let data: PrefetchData = if force || !data_file.exists() {
                            let sample_reference = plate
                                .get_sample_reference(&sample_id)
                                .ok_or_else(|| anyhow::anyhow!("No reference for {}", sample_id))?;

                            let (tags, ignore_taxstr, exclude_lod) = get_note_instructions(&sample_reference);

                            let tiered_filter_config = match ignore_taxstr {
                                Some(ignore_taxstr) => tiered_filter_config.with_ignore_taxstr(ignore_taxstr),
                                None => tiered_filter_config,
                            };

                            let config = MetaGpConfig::new(
                                sample_reference.sample_id.clone(),
                                sample_reference.sample_type.clone(),
                                sample_reference.result.clone(),
                                sample_reference.positive_taxa(),
                                exclude_lod,
                                Some(negative_controls),
                                Some(tags),
                                tiered_filter_config,
                                contam_config.clone(),
                            );

                            plate.prefetch(&client, &data_file, &config, prevalence_contamination)?
                        } else {
                            let fh = std::fs::File::open(&data_file)?;
                            serde_json::from_reader(fh)?
                        };

                        // Build JSON summary item
                        let counts = counts_by_category(&data);
                        let contam_counts = counts_by_category_contam(&data);
                        let (positive, positive_match) = positive_candidate_match(&data);

                        let summary = PerSampleSummary {
                            sample: data.config.sample.clone(),
                            counts,
                            contamination_counts: contam_counts,
                            positive,
                            positive_match,
                        };

                        // Build missed-detection row if applicable
                        let missed = if is_missed_detection(&data) {
                            Some(MissedDetectionRow {
                                sample: data.config.sample.clone(),
                                reference: reference_names_from_config(&data.config),
                                status: PrefetchStatus::NotDetected,
                            })
                        } else {
                            None
                        };

                        Ok(Some((summary, missed)))
                    })();

                    if let Err(e) = &result {
                        log::error!("Error prefetching {}: {:?}", sample_id, e);
                    }
                    result.unwrap_or(None)
                })
                .collect();

            // Split vectors
            let mut summaries: Vec<PerSampleSummary> = Vec::with_capacity(results.len());
            let mut missed_rows: Vec<MissedDetectionRow> = Vec::new();
            for (s, m) in results {
                summaries.push(s);
                if let Some(row) = m {
                    missed_rows.push(row);
                }
            }

            // Write JSON summary if requested
            if let Some(summary_path) = &args.summary {
                let total_positive = summaries.iter().filter(|s| s.positive).count();
                let total_positive_with_candidate_match = summaries.iter().filter(|s| s.positive_match).count();

                let overall = OverallSummary {
                    total_samples: summaries.len(),
                    total_positive,
                    total_positive_with_candidate_match,
                    per_sample: summaries,
                };

                if let Some(parent) = summary_path.parent() {
                    std::fs::create_dir_all(parent)?;
                }
                let mut fh = std::fs::File::create(summary_path)?;
                serde_json::to_writer_pretty(&mut fh, &overall)?;
                log::info!("Wrote summary JSON to {}", summary_path.display());
            }

            // Write TSV of missed detections if requested
            if let Some(missed_path) = &args.summary_missed {
                if let Some(parent) = missed_path.parent() {
                    std::fs::create_dir_all(parent)?;
                }
                write_tsv(&missed_rows, missed_path, true)?;
                log::info!("Wrote missed-detections TSV to {}", missed_path.display());
            }

            if !missed_rows.is_empty() {
                if let (Some(pd_table_path), Some(pd_out_path)) =
                    (&args.pathogen_detection_table, &args.pathogen_missed)
                {
            
                    let mut missed_map: HashMap<String, HashSet<String>> = HashMap::new();
                    for row in &missed_rows {
                        let refs = row
                            .reference
                            .split(';')
                            .map(|s| s.trim().to_string())
                            .filter(|s| !s.is_empty())
                            .collect::<HashSet<_>>();
                        missed_map
                            .entry(row.sample.clone())
                            .and_modify(|set| set.extend(refs.clone()))
                            .or_insert(refs);
                    }
            
                    // Read PD table TSV
                    let file = File::open(pd_table_path).map_err(|e| {
                        anyhow::anyhow!("Failed to open pathogen detection table {}: {e}", pd_table_path.display())
                    })?;
                    let mut rdr = csv::ReaderBuilder::new()
                        .delimiter(b',')
                        .has_headers(true)
                        .from_reader(BufReader::new(file));
            
                    let mut subset: Vec<PathogenDetectionTableRecord> = Vec::new();
            
                    for rec in rdr.deserialize::<PathogenDetectionTableRecord>() {
                        match rec {
                            Ok(row) => {
                                let lineage_str = match &row.lineage {
                                    Some(l) => l,
                                    None => continue,
                                };
                    
                                let mut keep = false;
                                for (sample, refs) in missed_map.iter() {
                                    if row.id.starts_with(sample)
                                        && refs.iter().any(|r| lineage_str.contains(r))
                                    {
                                        keep = true;
                                        break;
                                    }
                                }
                    
                                if keep {
                                    subset.push(row);
                                }
                            }
                            Err(err) => {
                                log::warn!("Skipping malformed PD table row: {err}");
                            }
                        }
                    }
            
                    // Write subset if requested
                    if !subset.is_empty() {
                        if let Some(parent) = pd_out_path.parent() {
                            std::fs::create_dir_all(parent)?;
                        }
                        write_tsv(&subset, pd_out_path, true)?;
                        log::info!(
                            "Wrote pathogen-missed subset (n={}) to {}",
                            subset.len(),
                            pd_out_path.display()
                        );
                    } else {
                        log::info!("No pathogen detection rows matched missed detections; no subset written.");
                    }
                }
            }
        }
        Commands::UploadCiqaDataset( args ) => {

            let api_client = CerebroClient::new(
                &cli.url,
                cli.token,
                false,
                cli.danger_invalid_certificate,
                cli.token_file,
                cli.team,
                cli.db,
                cli.project,
            )?;
            let fs_client = FileSystemClient::new(
                &api_client, 
                &cli.fs_url, 
                &cli.fs_port
            );

            for file in args.fastq_pe {
                let sample_id = get_file_component(&file, FileComponent::FileStem);
                fs_client.upload_files(
                    &Vec::from([file.to_path_buf()]), Some(
                        format!(
                            "ciqa-{}", 
                            args.run_id
                                .clone()
                                .unwrap_or(
                                    uuid::Uuid::new_v4().to_string()
                                )
                            )
                        ),
                        sample_id.ok(),
                        args.pipeline_id.clone(),
                        args.description.clone(), 
                        Some(FileType::ReadPaired),
                        UploadConfig::default(), 
                        None
                    )?;

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

                let mut reference_plate = ReferencePlate::from_review(
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
                    args.reference.clone(),
                    args.header_text.as_deref(),
                    Some(&consensus_stats)
                )?;
            }
        },

        #[cfg(feature = "local")]
        Commands::DiagnoseLocal( args ) => {

            std::fs::create_dir_all(&args.outdir)?;

            let state_dir = args.outdir.join("state_logs");
            std::fs::create_dir_all(&state_dir)?;

            let plate = ReferencePlate::from_review(
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

                let args = Arc::clone(&args);
                let plate = Arc::clone(&plate);
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
                            args.top_k,
                            args.top_p,
                            args.min_p,
                            gpu_id,
                        )
                    )?;

                    for sample_id in batch_samples {
                        
                        // rebuild paths

                        let data_file   = args.prefetch.join(format!("{sample_id}.prefetch.json"));

                        let result_file = args.outdir.join(format!("{sample_id}.model.json"));
                        let state_file  = state_dir.join(format!("{sample_id}.state.json"));
                        let bench_file  = state_dir.join(format!("{sample_id}.bench.json"));
                        
                        if !data_file.exists() {
                            log::info!("no prefetch for {}, skipping", sample_id);
                            continue;
                        }

                        if result_file.exists() && !args.force {
                            log::warn!("result file exists and force not activated - skipping {}", sample_id)
                        }
                        
                        // load prefetch
                        let prefetch_data = PrefetchData::from_json(&data_file)?;

                        // instantiate agent
                        let mut agent = DiagnosticAgent::new(args.task_config.clone(), args.tree_config.clone())?;
                        
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
                        
                        let post_filter = if args.post_filter { 
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
                            prefetch_data,
                            &mut text_generator,
                            if args.sample_context.unwrap_or(false) { Some(sample_context) } else { None },
                            if args.clinical_notes { clinical_notes } else { None },
                            args.assay_context.clone(),
                            args.agent_primer.clone(),
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

        Commands::McnemarAdjust( args ) => {

            let base = DiagnosticData::from_json(args.baseline)?;
            let base_reviews = base.get_consensus_reviews_filtered()
                .ok_or_else(|| anyhow::anyhow!("baseline missing consensus reviews"))?;
    
            // load ablations
            let mut labeled = Vec::new();
            for path in &args.comparisons {
                let dd = DiagnosticData::from_json(path)?;
                if let Some(revs) = dd.get_consensus_reviews_filtered() {
                    // use filename as ID
                    labeled.push((get_file_component(path, FileComponent::FileStem)?, revs));
                }
            }
    
            let rows = mcnemar_batch_adjust(&base_reviews, labeled, args.method);
    
            // print concise table
            println!("id\tb\tc\tn_discord\tp_raw\tp_adj\tmethod\tmode");
            for r in rows {
                println!("{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:?}\t{}",
                    r.id,
                    r.result.b,
                    r.result.c,
                    r.result.total_discordant,
                    r.p_raw,
                    r.p_adj,
                    args.method,
                    r.result.test_type
                );
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
        Commands::Tui( args ) => {

            start_tui()?;

        }
    }

    Ok(())

}


fn get_note_instructions(sample_reference: &SampleReference) -> (Vec<String>, Option<Vec<String>>, Option<bool>) {

    let tags = match sample_reference.note {
        Some(ref note) => {
            if note.contains("ignore_rna") {
                vec![String::from("DNA")]
            } else if note.contains("ignore_dna") {
                vec![String::from("RNA")]
            } else if note.contains("only_dna") {
                vec![String::from("DNA")]
            } else if note.contains("only_rna") {
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

    let exclude_lod = match sample_reference.note {
        Some(ref note) => Some(note.contains("exclude_lod")),
        None => None
    };

    (tags, ignore_taxstr, exclude_lod)
}
