#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unreachable_code)]


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

use cerebro_report::utils::init_logger;
use cerebro_report::report::ClinicalReport;
use cerebro_report::terminal::{App, Commands};

// NB: import error arises from including `cli.rs` in `lib.rs` for obvious reasons

fn main() -> anyhow::Result<()> {
    init_logger();

    let cli = App::parse();

    match &cli.command {
        Commands::Compile( args ) => {

            let clinical_report = ClinicalReport::from_toml(
                &args.base_config, 
                args.patient_config.clone(), 
                args.patient_header_config.clone(), 
                args.patient_result_config.clone(), 
                args.sample_config.clone()
            )?;

            if args.pdf {
                clinical_report.render_pdf(Some(args.output.clone()), None, None)?;
            } else {
                clinical_report.render_template(Some(args.output.clone()))?;
            }
            std::process::exit(0)
        }
    }
    Ok(())
}
