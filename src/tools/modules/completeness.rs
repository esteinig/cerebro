use anyhow::Result;
use std::io::BufWriter;
use std::{path::Path, fs::File};
use needletail::{parse_fastx_file, Sequence};
use std::str::from_utf8;
use crate::utils::CompressionExt;
use crate::tools::error::ToolError;

/*
========================
Custom error definitions
========================
*/


pub struct Completeness {
    output_format: Option<niffler::compression::Format>,
    compression_level: niffler::compression::Level
}

impl Completeness {
    pub fn new(
        output_format: Option<niffler::compression::Format>,
        compression_level: niffler::compression::Level
    ) -> Result<Self, ToolError> {
        Ok(Self {
            output_format, compression_level
        })
    }
    pub fn completeness(
        &self,
        input: &Path,
        output: &Path,
        min_completeness: &f64,
    ) -> Result<(), ToolError>{
        // Input output of read files includes compression detection
        let mut reader = parse_fastx_file(input)?;
        
       
        let file = File::create(&output)?;
        let file_handle = BufWriter::new(file);

        let fmt = match self.output_format {
            None => niffler::Format::from_path(&output),
            Some(format) => format,
        };
        let mut writer = niffler::get_writer(Box::new(file_handle), fmt, self.compression_level)?;

        while let Some(record) = reader.next() {
            let rec = record?;
            let count = rec.sequence().iter().filter(|x| *x == &b'N').count();
            let rec_id = from_utf8(rec.id())?.split(' ').next().unwrap_or("");
            let rec_len = rec.num_bases();
            let completeness  = 100.0 - (count as f64 / rec_len as f64)*100.0;

            if &completeness >= min_completeness {
                needletail::parser::write_fasta(rec.id(), rec.sequence(), &mut writer, needletail::parser::LineEnding::Unix)?;
            }
            
            println!("{:}\t{:}\t{:}\t{:.2}", rec_id, count, rec_len, completeness)
        }
        Ok(())
    }

}