use anyhow::Result;
use clap::Parser;

use cerebro_client::client::CerebroClient;
use cerebro_fs::capture::{CaptureRule, CaptureStatus};
use cerebro_fs::client::FileSystemClient;
use cerebro_fs::config::FsConfig;
use cerebro_fs::terminal::{App, Commands};
use cerebro_fs::{client::UploadConfig, utils::init_logger, weed::download_and_install_weed};

fn main() -> Result<()> {
    init_logger();

    let cli = App::parse();

    let api_client = CerebroClient::new(
        &cli.url,
        cli.token,
        match &cli.command {
            Commands::Login(_) | Commands::Ping(_) => true,
            _ => false,
        },
        cli.danger_invalid_certificate,
        cli.token_file,
        cli.team,
        cli.db,
        cli.project,
    )?;
    let fs_config = FsConfig {
        master_url: cli.fs_url.clone(),
        master_port: cli.fs_port.clone(),
        localhost: true,
        filer_url: cli.fs_filer_url.clone(),
        access: cli.fs_access.clone(),
        danger_invalid_certificate: cli.danger_invalid_certificate,
    };
    let fs_client = FileSystemClient::with_config(&api_client, fs_config);

    match &cli.command {
        Commands::Login(args) => {
            api_client.login_user(&args.email, args.password.clone(), false)?
        }
        Commands::Ping(_) => {
            fs_client.ping_status()?;
            api_client.ping_status()?;
        }
        Commands::Upload(args) => {
            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!(
                "Checking status of SeaweedFS master at {}",
                &fs_client.get_url()
            );
            fs_client.ping_status()?;

            log::info!("Starting file processing and upload");
            fs_client.upload_files(
                &args.files,
                args.run_id.clone(),
                args.sample_id.clone(),
                args.pipeline_id.clone(),
                args.description.clone(),
                args.file_type.clone(),
                UploadConfig {
                    tier: args.tier,
                    retention: args.retention,
                    legal_hold: args.legal_hold,
                    ..UploadConfig::default()
                },
                None,
            )?;
        }
        Commands::Download(args) => {
            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!(
                "Checking status of SeaweedFS master at {}",
                &fs_client.get_url()
            );
            fs_client.ping_status()?;

            log::info!("Starting file download");
            let report = fs_client.download_files(
                &args.fids,
                args.run_id.clone(),
                args.sample_id.clone(),
                &args.outdir,
                args.verify,
            )?;
            log::info!(
                "Downloaded {} file(s) to {}",
                report.written.len(),
                args.outdir.display()
            );
            if report.restore_required() {
                log::warn!(
                    "{} file(s) require an archival restore before download; run `cerebro-fs restore` first: {:?}",
                    report.restore_pending.len(),
                    report.restore_pending
                );
            }
        }
        Commands::Restore(args) => {
            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            let outcomes = fs_client.restore_files(args.run_id.clone(), args.sample_id.clone())?;
            let pending = outcomes
                .iter()
                .filter(|o| {
                    o.state == cerebro_model::api::files::retention::RestoreState::Requested
                })
                .count();
            for outcome in &outcomes {
                log::info!(
                    "{}: {} ({})",
                    outcome.identifier,
                    outcome.state,
                    outcome.message
                );
            }
            log::info!(
                "{} of {} file(s) require an archival restore",
                pending,
                outcomes.len()
            );
        }
        Commands::Capture(args) => {
            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!(
                "Checking status of SeaweedFS master at {}",
                &fs_client.get_url()
            );
            fs_client.ping_status()?;

            let upload_config = UploadConfig {
                retention_policy: cerebro_model::api::files::retention::RetentionPolicy::from_env(),
                ..UploadConfig::default()
            };
            let rules = CaptureRule::default_ruleset();

            log::info!(
                "Capturing pipeline outputs from {}",
                args.output_dir.display()
            );
            let report = fs_client.capture_outputs(
                args.run_id.clone(),
                args.sample_id.clone(),
                args.pipeline_id.clone(),
                &args.output_dir,
                &rules,
                &upload_config,
            )?;

            for outcome in &report.outcomes {
                match &outcome.status {
                    CaptureStatus::Captured { ftype, retention } => log::info!(
                        "captured {} [{:?}, {:?}]",
                        outcome.relative_path,
                        ftype,
                        retention
                    ),
                    CaptureStatus::Ignored => log::debug!("ignored {}", outcome.relative_path),
                    CaptureStatus::Failed(err) => {
                        log::warn!("FAILED {}: {}", outcome.relative_path, err)
                    }
                }
            }
            log::info!(
                "Capture complete: {} captured, {} ignored, {} failed",
                report.captured(),
                report.ignored(),
                report.failed()
            );
            if !report.ok() {
                log::warn!("Some pipeline outputs failed to capture");
            }
        }
        Commands::Manifest(args) => {
            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!(
                "Checking status of SeaweedFS master at {}",
                &fs_client.get_url()
            );
            fs_client.ping_status()?;

            // Provenance metadata: load from --metadata if given, then let the
            // explicit pipeline name/version flags override.
            let mut provenance = match &args.metadata {
                Some(path) => {
                    let bytes = std::fs::read(path).map_err(|e| {
                        anyhow::anyhow!("failed to read metadata file {}: {e}", path.display())
                    })?;
                    serde_json::from_slice::<cerebro_model::api::files::manifest::ManifestProvenance>(&bytes)
                        .map_err(|e| anyhow::anyhow!("invalid provenance metadata file: {e}"))?
                }
                None => cerebro_model::api::files::manifest::ManifestProvenance::default(),
            };
            if let Some(name) = &args.pipeline_name {
                provenance.pipeline_name = name.clone();
            }
            if let Some(version) = &args.pipeline_version {
                provenance.pipeline_version = version.clone();
            }

            let upload_config = UploadConfig {
                retention_policy: cerebro_model::api::files::retention::RetentionPolicy::from_env(),
                ..UploadConfig::default()
            };

            log::info!("Building run provenance manifest");
            let manifest = fs_client.build_run_manifest(
                args.run_id.clone(),
                args.sample_id.clone(),
                args.pipeline_id.clone(),
                provenance,
            )?;
            log::info!(
                "Manifest for pipeline {} {} : {} input(s), {} output(s), seal {}",
                manifest.pipeline_name,
                manifest.pipeline_version,
                manifest.inputs.len(),
                manifest.outputs.len(),
                manifest.content_hash.as_deref().unwrap_or("none")
            );

            let id = fs_client.capture_manifest(&manifest, &upload_config)?;
            log::info!("Captured run manifest as file {}", id);
        }
        Commands::ReportOut(args) => {
            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            // Anchor timestamp: explicit RFC 3339, else now.
            let reported_at = match &args.reported_at {
                Some(ts) => chrono::DateTime::parse_from_rfc3339(ts)
                    .map_err(|e| anyhow::anyhow!("invalid --reported-at timestamp '{ts}': {e}"))?
                    .with_timezone(&chrono::Utc),
                None => chrono::Utc::now(),
            };

            // Retention durations are deployment configuration; the default policy
            // is used for this preview.
            let policy = cerebro_model::api::files::retention::RetentionPolicy::from_env();

            log::info!(
                "Planning report-out lifecycle at {}",
                reported_at.to_rfc3339()
            );
            let report = fs_client.plan_report_out(
                args.run_id.clone(),
                args.sample_id.clone(),
                reported_at,
                &policy,
                args.warm,
            )?;
            for entry in &report.entries {
                let until = entry
                    .transition
                    .retain_until
                    .map(|t| t.to_rfc3339())
                    .unwrap_or_else(|| "indefinite".to_string());
                match entry.transition.cold_move_at {
                    Some(cold_at) => log::info!(
                        "{}: -> tier {} (warm→cold at {}) , retain until {}",
                        entry.name,
                        entry.transition.target_tier,
                        cold_at.to_rfc3339(),
                        until
                    ),
                    None => log::info!(
                        "{}: -> tier {} , retain until {}",
                        entry.name,
                        entry.transition.target_tier,
                        until
                    ),
                }
            }

            if args.persist {
                use cerebro_model::api::files::schema::UpdateFileLifecycleSchema;
                log::info!(
                    "Persisting report-out lifecycle for {} file(s)",
                    report.entries.len()
                );
                for entry in &report.entries {
                    let mut schema = UpdateFileLifecycleSchema {
                        tier: Some(entry.transition.target_tier),
                        reported_at: Some(entry.transition.reported_at),
                        ..Default::default()
                    };
                    match entry.transition.retain_until {
                        Some(retain_until) => schema.retain_until = Some(retain_until),
                        None => schema.clear_retain_until = true,
                    }
                    api_client.update_file_lifecycle(&entry.id, &schema)?;
                }
                log::info!("Report-out lifecycle persisted (physical tier moves are applied by the worker)");
            } else {
                log::info!(
                    "Planned report-out for {} file(s) (preview; re-run with --persist to write, tier moves are applied by the worker)",
                    report.entries.len()
                );
            }
        }
        Commands::Verify(args) => {
            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!(
                "Checking status of SeaweedFS master at {}",
                &fs_client.get_url()
            );
            fs_client.ping_status()?;

            log::info!(
                "Verifying file integrity{}",
                if args.repair { " (repair enabled)" } else { "" }
            );
            let report =
                fs_client.verify_files(args.run_id.clone(), args.sample_id.clone(), args.repair)?;
            for outcome in &report.outcomes {
                let replicas = outcome
                    .replicas
                    .map(|n| n.to_string())
                    .unwrap_or_else(|| "?".to_string());
                log::info!(
                    "{} [{} replica(s)]: {:?}",
                    outcome.name,
                    replicas,
                    outcome.status
                );
            }
            log::info!(
                "Integrity: {} ok, {} failed of {} file(s)",
                report.ok_count(),
                report.failed_count(),
                report.outcomes.len()
            );
            if !report.ok() {
                log::warn!("Integrity verification found unrecoverable mismatches");
            }
        }
        Commands::Health => {
            let health = fs_client.topology_health();
            for component in &health.components {
                if component.reachable {
                    log::info!("{}: ok ({})", component.component, component.detail);
                } else {
                    log::warn!(
                        "{}: UNREACHABLE ({})",
                        component.component,
                        component.detail
                    );
                }
            }
            log::info!(
                "Cerebro FS topology health: {}",
                if health.healthy() {
                    "healthy"
                } else {
                    "degraded"
                }
            );
        }
        Commands::Stage(args) => {
            fs_client.stage_files(&args.json, &args.outdir, args.pipeline.clone())?;
        }
        Commands::Hash(args) => {
            let n = cerebro_fs::hash::hash_directory(&args.outdir, &args.output)?;
            log::info!(
                "Wrote {} BLAKE3 checksum(s) to {}",
                n,
                args.output.display()
            );
        }
        Commands::Delete(args) => {
            fs_client.delete_files(
                &args.file_ids,
                args.run_id.clone(),
                args.sample_id.clone(),
                args.all,
            )?;
        }
        Commands::List(args) => {
            api_client.list_files(
                args.run_id.clone(),
                args.watcher_id.clone(),
                args.page,
                args.limit,
                true,
            )?;
        }
        Commands::Weed(args) => {
            download_and_install_weed(&args.version, &args.outdir.join("weed"))?
        }
    }

    Ok(())
}
