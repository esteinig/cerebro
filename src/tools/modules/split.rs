use anyhow::Result;
use std::fs::File;
use needletail::parser::{write_fasta, LineEnding};
use std::path::{Path, PathBuf};
use needletail::{parse_fastx_file, Sequence};
use std::io::BufWriter;
use crate::utils::CompressionExt;
use crate::tools::error::ToolError;


pub struct Splitter {
    pub outdir: PathBuf,
    pub output_format: Option<niffler::compression::Format>,
    pub compression_level: niffler::compression::Level
}

impl Splitter {
    pub fn new(
        outdir: &PathBuf,
        output_format: &Option<niffler::compression::Format>,
        compression_level: &niffler::compression::Level
    ) -> Result<Self, ToolError> {
        std::fs::create_dir_all(&outdir)?;
        Ok(Self { outdir: outdir.to_owned(), output_format: output_format.to_owned(), compression_level: compression_level.to_owned() })
    }
    pub fn split(
        &self,
        input: &Path,
    ) -> Result<(), ToolError>{
        // Input output of read files includes compression detection
        let mut reader = parse_fastx_file(input)?;

        // Automatically uses input file extension
        let file_name = match input.file_name() {
            Some(fname) => match fname.to_str(){
                Some(_fname) => _fname,
                None => return Err(ToolError::FileName)
            },
            None => return Err(ToolError::FileName)
        };
        
        let mut idx = 0;
        while let Some(record) = reader.next() {
            let rec = record?;
            let mut rec_path = self.outdir.clone();
            rec_path.push(format!("{:0>2}-{:}", idx, file_name));

            let file = File::create(&rec_path)?;
            let file_handle = BufWriter::new(file);
            let fmt = match self.output_format {
                None => niffler::Format::from_path(&rec_path),
                Some(format) => format,
            };
            let mut writer = niffler::get_writer(Box::new(file_handle), fmt, self.compression_level)?;
            
            write_fasta(rec.id(), rec.sequence(), &mut writer, LineEnding::Unix)?;  // only unix
            idx += 1;
        }
        Ok(())
    }

}