// Cerebro: metagenomic and -transcriptomic diagnostics for clinical production environments
// Copyright (C) 2023  Eike Steinig

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

use clap::Parser;
use crate::utils::init_logger;
use crate::terminal::{App, Commands};
use crate::slack::{SlackMessage, SlackMessenger};
use crate::watcher::watch_production;
use crate::watcher::SlackConfig;

mod error;
mod slack;
mod terminal;
mod launcher;
mod watcher;
mod utils;


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

            if let Err(error) = watch_production(
                &args.path, 
                std::time::Duration::from_secs(args.interval),
                std::time::Duration::from_secs(args.timeout),
                std::time::Duration::from_secs(args.timeout_interval),
                SlackConfig { channel: args.slack_channel.clone(), token: args.slack_token.clone() }
            ) {
                log::error!("Error: {error:?}");
            }
        },
    }

    Ok(())

}
   