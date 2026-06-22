use anyhow::Result;
use cerebro_client::client::CerebroClient;
use cerebro_fs::client::FileSystemClient;
use cerebro_fs::config::FsConfig;
use cerebro_tower::{
    client::TowerClient,
    terminal::{App, Commands},
    tower::{CerebroTower, NextflowConfig},
    utils::init_logger,
};
use clap::Parser;
use tokio::main;

#[tokio::main]
async fn main() -> Result<()> {
    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Slack(args) => {
            // Handle Slack command (implement your logic here)
        }
        Commands::Watch(args) => {
            let client = TowerClient::new(
                &cli.url,
                cli.token.clone(),
                cli.danger_invalid_certificate,
                cli.token_file.clone(),
                cli.team.clone(),
                cli.db.clone(),
                cli.project.clone(),
            )?;

            client.ping_servers().await?;

            // Cerebro FS client: the tower constructs a FileSystemClient
            // so it can stage inputs and — in Stage 2 — capture and upload
            // pipeline output artefacts. The same API credentials are reused.
            let api_client = CerebroClient::new(
                &cli.url,
                cli.token.clone(),
                false,
                cli.danger_invalid_certificate,
                cli.token_file.clone(),
                cli.team.clone(),
                cli.db.clone(),
                cli.project.clone(),
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
            // ping_status() is blocking (reqwest::blocking); run it off the
            // async runtime to avoid nested-runtime panics.
            let fs_for_ping = fs_client.clone();
            tokio::task::spawn_blocking(move || fs_for_ping.ping_status()).await??;

            let nextflow = NextflowConfig::new(
                &args.main,
                &args.config,
                &args.workdir,
                &args.databases,
                &args.profile,
                args.cleanup,
            )
            .await?;

            let tower = CerebroTower::new(client, nextflow, fs_client)?;

            tower.watch(&args.tower_id, true).await?;
        }
    }

    Ok(())
}
