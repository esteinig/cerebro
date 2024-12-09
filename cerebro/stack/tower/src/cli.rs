

use cerebro_tower::{client::TowerClient, terminal::{App, Commands}, tower::{CerebroTower, NextflowConfig}, utils::init_logger};
use clap::Parser;
use anyhow::Result;
use tokio::main;

#[tokio::main]
async fn main() -> Result<()> {

    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Slack(args) => {
            // Handle Slack command (implement your logic here)
        },
        Commands::Watch(args) => {


            let client = TowerClient::new(
                &cli.url,
                cli.token,
                cli.danger_invalid_certificate,
                cli.token_file,
                cli.team,
                cli.db,
                cli.project,
            )?;

            client.ping_servers().await?;

            let nextflow = NextflowConfig::new(
                &args.main, 
                &args.config, 
                &args.workdir, 
                &args.databases,
                &args.profile,
                args.cleanup
            ).await?;

            let tower = CerebroTower::new(
                client, 
                nextflow
            )?;

            tower.watch(&args.tower_id, true).await?;
        },
    }

    Ok(())
}