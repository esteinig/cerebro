
use clap::Parser;
use cerebro_watcher::utils::init_logger;
use cerebro_watcher::terminal::{App, Commands};
use cerebro_watcher::slack::{SlackMessage, SlackMessenger};
use cerebro_watcher::watcher::watch_production;
use cerebro_watcher::watcher::SlackConfig;



fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Slack( args ) => {
            let messenger = SlackMessenger::new(&args.token);
            messenger.send(
                &SlackMessage::new(&args.channel, &args.message)
            ).expect("Failed to send message");
        },
        Commands::Watch( args ) => {

            let slack_config = match (&args.slack_channel, &args.slack_token) {
                (Some(channel), Some(token)) => Some(
                    SlackConfig { channel: channel.to_string(), token: token.to_string() }
                ),
                _ => None
            };

            if let Err(error) = watch_production(
                &args.path, 
                std::time::Duration::from_secs(args.interval),
                std::time::Duration::from_secs(args.timeout),
                std::time::Duration::from_secs(args.timeout_interval),
                slack_config
            ) {
                log::error!("Error: {error:?}");
            }
        },
    }

    Ok(())

}
   