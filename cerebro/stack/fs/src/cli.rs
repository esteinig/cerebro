
use clap::Parser;
use anyhow::Result;

use cerebro_fs::{client::UploadConfig, utils::init_logger};
use cerebro_model::api::files::model::WatcherConfig;
use cerebro_fs::client::FileSystemClient;
use cerebro_client::client::CerebroClient;
use cerebro_fs::terminal::{App, Commands};


fn main() -> Result<()> {

    init_logger();

    let cli = App::parse();

    let api_client = CerebroClient::new(
        &cli.url,
        &cli.token,
        match &cli.command { Commands::Login(_) | Commands::Ping(_) => true, _ => false },
        cli.danger_invalid_certificate,
        &cli.token_file
    )?;
    let fs_client = FileSystemClient::new(
        &api_client, 
        &cli.fs_url, 
        &cli.fs_port
    );

    match &cli.command {
        Commands::Login( args ) => {
            api_client.login_user(&args.email, &args.password)?
        },
        Commands::Ping(_) => {
            fs_client.ping_status()?;
            api_client.ping_status()?;
        },           
        Commands::Upload( args ) => {

            log::info!("Checking status of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!("Checking status of SeaweedFS master at {}", &cli.fs_url);
            fs_client.ping_status()?;

            log::info!("Starting file processing and upload");
            fs_client.upload_files(
                &args.files,  
                &args.team_name,
                &args.db_name,
                &args.run_id,
                &args.sample_id,
                UploadConfig::default(),  // TODO
                WatcherConfig::default()  // TODO
            )?;
        },
        Commands::Download( args ) => {
            
        },
        Commands::Stage( args ) => {
            
        },
        Commands::Delete( args ) => {
            fs_client.delete_files(
                &args.file_ids,  
                &args.team_name,
                &args.db_name,
                args.run_id.clone()
            )?;
        },
        Commands::List( args ) => {
            api_client.get_files(
                &args.team_name,
                &args.db_name,
                args.run_id.clone(),
                args.page,
                args.limit,
                true
            )?;
        },

    }

    Ok(())
}