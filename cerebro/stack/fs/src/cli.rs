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
        match &cli.command { Commands::Login(_) | Commands::PingApi(_) => true, _ => false },
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
        Commands::PingApi(_) => {
            api_client.ping_status()?
        },           
        Commands::Upload( args ) => {

            log::info!("Checking availability of authenticated Cerebro API {}", &cli.url);
            api_client.ping_servers()?;

            log::info!("Checking availability of SeaweedFS master at {}", &cli.fs_url);
            fs_client.ping_status()?;

            log::info!("Starting file processing and upload");
            fs_client.upload_files(
                &args.files,  
                &args.team_name,
                &args.db_name,
                UploadConfig::default(),
                WatcherConfig::default()
            )?;
        },
        Commands::Download( args ) => {
            
        },
        Commands::List( args ) => {
            api_client.list_files(
                &args.team_name,
                &args.db_name,
                args.page,
                args.limit,
            )?;
        },

    }

    Ok(())
}