
use std::time::Duration;

use cerebro_client::client::CerebroClient;
use cerebro_fs::client::{FileSystemClient, UploadConfig};
use cerebro_model::api::watchers::model::ProductionWatcher;
use cerebro_watcher::utils::WatcherConfigArgs;
use cerebro_watcher::watcher::CerebroWatcher;
use cerebro_watcher::utils::{init_logger, UploadConfigArgs};
use cerebro_watcher::terminal::{App, Commands};
use cerebro_watcher::slack::{SlackMessage, SlackClient};
use cerebro_watcher::slack::SlackConfig;
use clap::Parser;

fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Slack( args ) => {
            let messenger = SlackClient::new(&args.token);
            messenger.send(
                &SlackMessage::new(&args.channel, &args.message)
            )?;
        },
        Commands::Watch( args ) => {


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

            let fs_client = FileSystemClient::new(
                &api_client, 
                &cli.fs_url, 
                &cli.fs_port
            );

            log::info!("Checking status of Cerebro API at {}",  &api_client.url);
            api_client.ping_servers()?;
    
            log::info!("Checking status of SeaweedFS master at {}",  &fs_client.fs_url);
            fs_client.ping_status()?;

            let watcher = CerebroWatcher::new(
                ProductionWatcher::from_args(args, &api_client)?, 
                api_client, 
                fs_client, 
                UploadConfig::from_args(args), 
                SlackConfig::from_args(args)
            )?;

            watcher.watch(
                &args.path, 
                Duration::from_secs(args.interval),
                Duration::from_secs(args.timeout),
                Duration::from_secs(args.timeout_interval),
            )?;

        },
    }

    Ok(())

}
   