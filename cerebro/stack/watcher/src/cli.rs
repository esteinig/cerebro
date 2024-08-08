
use cerebro_fs::client::UploadConfig;
use cerebro_model::api::files::model::WatcherConfig;
use cerebro_watcher::watcher::{CerebroClientConfig, CerebroWatcher};
use cerebro_watcher::utils::{init_logger, UploadConfigArgs};
use cerebro_watcher::terminal::{App, Commands};
use cerebro_watcher::slack::{SlackMessage, SlackClient};
use cerebro_watcher::slack::SlackConfig;
use cerebro_watcher::utils::{WatcherConfigArgs, CerebroClientConfigArgs};
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

            let slack_config = match (&args.slack_channel, &args.slack_token) {
                (Some(channel), Some(token)) => Some(
                    SlackConfig { channel: channel.to_string(), token: token.to_string() }
                ),
                _ => {
                    log::info!("No slack configuration provided to watcher");
                    None
                }
            };

            let client_config = CerebroClientConfig::from_args(&cli);
            let watcher_config = WatcherConfig::from_args(args);
            let upload_config = UploadConfig::from_args(args);

            let watcher = CerebroWatcher::new(watcher_config, client_config, upload_config, slack_config)?;
            watcher.watch(&args.path, args.glob.clone())?;

        },
    }

    Ok(())

}
   