
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
            let written = fs_client.download_files(
                &args.fids,
                args.run_id.clone(),
                args.sample_id.clone(),
                &args.outdir,
                args.verify,
            )?;
            log::info!("Downloaded {} file(s) to {}", written.len(), args.outdir.display());
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