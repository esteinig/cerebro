
use std::{fs::create_dir_all, path::PathBuf};

use cerebro_gp::{error::GptError, gpt::{AssayContext, SampleContext, DiagnosticAgent, DiagnosticResult, PrefetchData}, text::{GeneratorConfig, TextGenerator}};
use cerebro_model::api::cerebro::schema::{CerebroIdentifierSchema, MetaGpConfig, PrevalenceOutliers};
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

    match &cli.command {
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
                let config_file = args.outdir.join(format!("{sample_id}.config.json"));

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
                        args.contam_history,
                        None
                    )?;

                    log::info!("{:#?}", gp_config);

                    DiagnosticAgent::new(Some(api_client.clone()))
                        .await?
                        .prefetch(&data_file, &gp_config)
                        .await?;
                    
                    gp_config.to_json(&config_file)?;

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
                    args,
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
                args,
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

            create_dir_all(&args.outdir)?;

            let plate = ReferencePlate::new(
                &args.plate, 
                None,
                MissingOrthogonal::Indeterminate,
                false
            )?;

            for sample_id in &plate.samples {

                let data_file = args.prefetch.join(format!("{sample_id}.prefetch.json"));

                let result_file = args.outdir.join(format!("{sample_id}.result.json"));
                let state_file = args.outdir.join(format!("{sample_id}.state.json"));

                if data_file.exists() {
                    let prefetch_data = PrefetchData::from_json(data_file)?;

                    log::info!("{:#?}", prefetch_data.config);

                    let mut agent = DiagnosticAgent::new(None).await?;

                    let mut text_generator = TextGenerator::new(
                        GeneratorConfig::with_default(
                            args.model.clone(),
                            args.model_dir.clone(),
                            args.sample_len,
                            args.temperature,
                            args.gpu
                        )
                    )?;
                    
                    let sample_reference = plate.get_sample_reference(sample_id);

                    let sample_context = match sample_reference {
                        Some(ref sample) => {
                            if sample.sample_type == SampleType::Csf {
                                SampleContext::Csf
                            } else if sample.sample_type == SampleType::Eye {
                                SampleContext::Eye
                            } else {
                                SampleContext::None
                            }
                        },
                        None => SampleContext::None
                    };

                    let clinical_notes = match sample_reference {
                        None => None,
                        Some(ref sample) => sample.clinical.clone()
                    };

                    let result = agent.run_local(
                        &mut text_generator, 
                        match args.sample_context { 
                            Some(include) => if include { sample_context } else { SampleContext::None },
                            None => SampleContext::None
                        }, 
                        if args.clinical_notes { clinical_notes } else { None }, 
                        args.assay_context.clone(), 
                        &prefetch_data.config.clone(), 
                        Some(prefetch_data)
                    ).await?;

                    result.to_json(&result_file)?;
                    log::info!("{:#?}", result);
                    agent.state.to_json(&state_file)?;

                } else {
                    log::info!("Data file '{}' does not exists skipping sample '{}'", data_file.display(), sample_id);
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