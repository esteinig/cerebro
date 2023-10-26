
use anyhow::Result;
use needletail::{parse_fastx_file, Sequence};
use needletail::parser::{write_fastq, LineEnding};
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::str::from_utf8;
use uuid::Uuid;
use std::collections::HashMap;
use std::fmt;
use crate::utils::CompressionExt;
use crate::tools::error::ToolError;

/*
========================
Fake Illumina header
========================
*/

#[derive(Debug, Clone)]
pub struct FakeIlluminaHeader {
    instrument_id: String,
    run_number: u32,
    flowcell_id: String,
    lane: u32,
    tile: u32,
    x_pos: u32,
    y_pos: u32,
    umi: String,
    read_number: u8,
    is_filtered: String,
    control_number: u8,
    index: String
}

impl FakeIlluminaHeader {
    pub fn new(read_number: u8, flowcell_id: String) -> Self {

        FakeIlluminaHeader {
            instrument_id: "FCKZCK".to_string(),
            run_number: 1,
            flowcell_id,
            lane: 69,
            tile: 69,
            x_pos: 6969,
            y_pos: 6969,
            umi: "GATTACT+GTCTTAAC".to_string(),
            read_number,
            is_filtered: "N".to_string(),
            control_number: 0,
            index: "ATCCGA".to_string() 
        }
    }
}

impl fmt::Display for FakeIlluminaHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f, "{}:{}:{}:{}:{}:{}:{}:{} {}:{}:{}:{}", 
            self.instrument_id,
            self.run_number,
            self.flowcell_id,
            self.lane,
            self.tile,
            self.x_pos,
            self.y_pos,
            self.umi,
            self.read_number,
            self.is_filtered,
            self.control_number,
            self.index
        )
    }
}


pub struct ReadAnonymizer {
    output_format: Option<niffler::compression::Format>,
    compression_level: niffler::compression::Level
}

impl ReadAnonymizer {
    pub fn new(
        output_format: Option<niffler::compression::Format>,
        compression_level: niffler::compression::Level
    ) -> Result<Self, ToolError> {
                
        Ok(Self { output_format, compression_level })
    }

    pub fn anonymize_single_end(
        &self, 
        input: &Path,
        output: &Path,
    ) -> Result<(), ToolError> {

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
            let new_id = Uuid::new_v4().to_string();
            
            write_fastq(new_id.as_bytes(), rec.sequence(), rec.qual(), &mut writer, LineEnding::Unix)?;  // only unix
        }
        Ok(())
    }
    pub fn anonymize_paired_end(
        &self, 
        input: &Vec<PathBuf>,
        output: &Vec<PathBuf>,
        fake_illumina_header: &bool
    ) -> Result<(), ToolError> {

        let mut read_uuids: HashMap<String, String> = HashMap::new();
        for (i, input_file) in input.iter().enumerate() {
            // Input output of read files includes compression detection
            let mut reader = parse_fastx_file(input_file)?;

            let file = File::create(&output[i])?;
            let file_handle = BufWriter::new(file);
            let fmt = match self.output_format {
                None => niffler::Format::from_path(&output[i]),
                Some(format) => format,
            };

            let mut writer = niffler::get_writer(Box::new(file_handle), fmt, self.compression_level)?;

            while let Some(record) = reader.next() {
                let rec = record?;
                let rec_id = from_utf8(rec.id())?.split(' ').next().unwrap_or(""); // needletail parses the entire header as identifier (including description), this is all terribly inefficient
                
                let new_id = if !fake_illumina_header {
                    // Implement UUID header logic
                    if i == 0 {
                        let _new_id = Uuid::new_v4().to_string();
                        read_uuids.insert(rec_id.to_string(), _new_id.clone());
                        _new_id
                    } else {
                        // Create new UUID if the record identifier in the second file is not present in the first
                        read_uuids.get(rec_id).unwrap_or(&Uuid::new_v4().to_string()).to_string()
                    }    
                } else {
                    // Implement fake Illumina header
                    if i == 0 {
                        let _new_id = Uuid::new_v4().simple().to_string();
                        read_uuids.insert(rec_id.to_string(), _new_id.clone());
                        let _new_header = FakeIlluminaHeader::new(1, _new_id);
                        format!("{}", _new_header)
                    } else {
                        // Create new UUID if the record identifier in the second file is not present in the first
                        let _uuid = read_uuids.get(rec_id).unwrap_or(&Uuid::new_v4().simple().to_string()).to_string();
                        let _new_header = FakeIlluminaHeader::new(2, _uuid);
                        format!("{}", _new_header)
                    }    
                };
                
                write_fastq(new_id.as_bytes(), rec.sequence(), rec.qual(), &mut writer, LineEnding::Unix)?;  // only unix
            }
        };
        
        Ok(())
    }

}
