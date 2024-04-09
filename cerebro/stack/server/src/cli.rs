
use clap::Parser;
use crate::utils::init_logger;
use crate::terminal::{App, Commands};

mod api;
mod utils;
mod terminal;

fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Run( _ ) => { api::server::main()? },
    }

    Ok(())

}
   