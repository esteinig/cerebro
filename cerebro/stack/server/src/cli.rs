use crate::terminal::{App, Commands};
use crate::utils::init_logger;
use clap::Parser;

mod api;
mod terminal;
mod utils;

fn main() -> anyhow::Result<()> {
    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Run(_) => api::server::main()?,
    }

    Ok(())
}
