// Cerebro: metagenomic and -transcriptomic diagnostics for clinical production environments
// Copyright (C) 2024  Eike Steinig

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

#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(unreachable_code)]
#![allow(unreachable_patterns)]

use std::fs::create_dir_all;
use std::time::Duration;

use cerebro_client::error::HttpClientError;
use cerebro_model::api::towers::schema::RegisterTowerSchema;
use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;
use cerebro_model::api::watchers::schema::RegisterWatcherSchema;

use cerebro_watcher::utils::WatcherConfigArgs;
use cerebro_watcher::utils::UploadConfigArgs;

use cerebro_pipeline::taxon::TaxonThresholdConfig;
use clap::Parser;
use rayon::prelude::*;

use cerebro::tools::password::hash_password;
use cerebro::utils::init_logger;

use cerebro::terminal::{App, Commands, StackCommands};


fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    match &cli.command {
        // Subcommands for stack configuration and deployment
        Commands::Stack( subcommand ) => {
            match subcommand {
                StackCommands::Deploy( args ) => {
                    
                    let mut stack = cerebro::stack::deploy::StackConfig::from_toml(&args.config)?;

                    stack.configure( 
                        &args.outdir, 
                        args.dev, 
                        args.subdomain.clone(), 
                        args.trigger, 
                        args.fs_primary.clone(), 
                        args.fs_secondary.clone()
                    )?;

                    stack.clone_and_checkout_repository_process(
                        &args.git_url, 
                        args.branch.clone(),
                        args.revision.clone()
                    )?;
                },
                StackCommands::HashPassword( args ) => {
                    println!("{}", hash_password(&args.password)?);
                },  
            }
        },
        _ => {}
    }

    Ok(())
}