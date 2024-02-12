use itertools::Itertools;
use rayon::prelude::*;
use std::{path::PathBuf, collections::{HashSet, HashMap}, str::from_utf8};
use needletail::{parse_fastx_file, parser::write_fastq};
use crate::{tools::{error::ToolError, utils::{get_compression_writer, get_seq_record_identifier}}, pipeline::utils::get_file_stem};


pub struct Crosstalk { 
    output_format: Option<niffler::compression::Format>,
    compression_level: Option<niffler::compression::Level>
}
impl Crosstalk {
    pub fn new(output_format: Option<niffler::compression::Format>, compression_level: Option<niffler::compression::Level>) -> Self {
        Self {
            output_format,
            compression_level
        }
    }
    
}