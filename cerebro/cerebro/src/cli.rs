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

// #![allow(dead_code)]
// #![allow(unused_imports)]
// #![allow(unused_variables)]
// #![allow(unreachable_code)]
// #![allow(unreachable_patterns)]

use clap::Parser;

use utils::init_logger;
use api::report::ClinicalReport;
use terminal::ReportCommands;

use terminal::{
    App, 
    Commands, 
    StackCommands, 
};

mod tools;
mod utils;
mod stack;
mod terminal;

fn main() -> anyhow::Result<()> {

    init_logger();

    let cli = App::parse();

    match &cli.command {
        
        Commands::Report( subcommand ) => {
            match subcommand {
                ReportCommands::Compile( args ) => {
                    let clinical_report = ClinicalReport::from_toml(
                        &args.base_config, 
                        args.patient_config.clone(), 
                        args.patient_header_config.clone(), 
                        args.patient_result_config.clone(), 
                        args.sample_config.clone()
                    )?;

                    #[cfg(feature = "pdf")]
                    {
                        if args.pdf {
                            clinical_report.render_pdf(Some(args.output.clone()), None)?;
                        } else {
                            clinical_report.render_template(Some(args.output.clone()))?;
                        }
                        std::process::exit(0)
                    }
                    clinical_report.render_template(Some(args.output.clone()))?;
                },
                _
            }
        },
        // Subcommands for stack configuration and deployment
        Commands::Stack( subcommand ) => {

            match subcommand {
                // Configure the stack deployment by creating the 
                // output directory with templated assets
                StackCommands::Deploy( args ) => {
                    let mut stack = stack::deploy::Stack::from_toml(&args.config)?;
                    stack.configure( &args.outdir, args.dev, args.subdomain.clone(), args.trigger )?;

                    #[cfg(feature = "libgit")]
                    {
                        stack.clone_and_checkout_repository_libgit(
                            &args.git_url, 
                            args.libgit,
                            args.branch.clone(),
                            args.revision.clone(), 
                            args.git_ssh_key.clone(), 
                            args.git_ssh_pwd.clone()
                        )?;
                        std::process::exit(0)  // premature exit ok

                    }
                    stack.clone_and_checkout_repository_process(
                        &args.git_url, 
                        args.branch.clone(),
                        args.revision.clone()
                    )?;
                    
                }
            }
        },
    }
}