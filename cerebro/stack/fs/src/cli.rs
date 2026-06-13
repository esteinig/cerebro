
use clap::Parser;
use anyhow::Result;

use cerebro_fs::{client::UploadConfig, utils::init_logger, weed::download_and_install_weed};
use cerebro_fs::client::FileSystemClient;
use cerebro_fs::config::FsConfig;
use cerebro_client::client::CerebroClient;
use cerebro_fs::terminal::{App, Commands};


fn main() -> Result<()> {

    init_logger();

    let cli = App::parse();

    let api_client = CerebroClient::new(
        &cli.url,
        cli.token,
        match &cli.command { Commands::Login(_) | Commands::Ping(_) => true, _ => false },
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
        Commands::Login( args ) => {
            api_client.login_user(&args.email, args.password.clone(), false)?
        },
        Commands::Ping(_) => {
            fs_client.ping_status()?;
            api_client.ping_status()?;
        },           
        Commands::Upload( args ) => {

            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!("Checking status of SeaweedFS master at {}", &fs_client.get_url());
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
        },
        Commands::Download( args ) => {

            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!("Checking status of SeaweedFS master at {}", &fs_client.get_url());
            fs_client.ping_status()?;

            log::info!("Starting file download");
            let report = fs_client.download_files(
                &args.fids,
                args.run_id.clone(),
                args.sample_id.clone(),
                &args.outdir,
                args.verify,
            )?;
            log::info!("Downloaded {} file(s) to {}", report.written.len(), args.outdir.display());
            if report.restore_required() {
                log::warn!(
                    "{} file(s) require an archival restore before download; run `cerebro-fs restore` first: {:?}",
                    report.restore_pending.len(),
                    report.restore_pending
                );
            }
        },
        Commands::Restore( args ) => {

            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            let outcomes = fs_client.restore_files(
                args.run_id.clone(),
                args.sample_id.clone(),
            )?;
            let pending = outcomes.iter().filter(|o| o.state == cerebro_model::api::files::retention::RestoreState::Pending).count();
            for outcome in &outcomes {
                log::info!("{}: {} ({})", outcome.identifier, outcome.state, outcome.message);
            }
            log::info!("{} of {} file(s) require an archival restore", pending, outcomes.len());
        },
        Commands::Verify( args ) => {

            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!("Checking status of SeaweedFS master at {}", &fs_client.get_url());
            fs_client.ping_status()?;

            log::info!("Verifying file integrity{}", if args.repair { " (repair enabled)" } else { "" });
            let report = fs_client.verify_files(
                args.run_id.clone(),
                args.sample_id.clone(),
                args.repair,
            )?;
            for outcome in &report.outcomes {
                let replicas = outcome.replicas.map(|n| n.to_string()).unwrap_or_else(|| "?".to_string());
                log::info!("{} [{} replica(s)]: {:?}", outcome.name, replicas, outcome.status);
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
        },
        Commands::Health => {

            let health = fs_client.topology_health();
            for component in &health.components {
                if component.reachable {
                    log::info!("{}: ok ({})", component.component, component.detail);
                } else {
                    log::warn!("{}: UNREACHABLE ({})", component.component, component.detail);
                }
            }
            log::info!("Cerebro FS topology health: {}", if health.healthy() { "healthy" } else { "degraded" });
        },
        Commands::Stage( args ) => {
            fs_client.stage_files(
                &args.json,
                &args.outdir,
                args.pipeline.clone()
            )?;
            
        },
        Commands::Delete( args ) => {
            fs_client.delete_files(
                &args.file_ids,  
                args.run_id.clone(),
                args.sample_id.clone(),
                args.all
            )?;
        },
        Commands::List( args ) => {
            api_client.list_files(
                args.run_id.clone(),
                args.watcher_id.clone(),
                args.page,
                args.limit,
                true
            )?;
        },
        Commands::Weed( args ) => {
            download_and_install_weed(
                &args.version, 
                &args.outdir.join("weed")
            )?
        },

    }

    Ok(())
}