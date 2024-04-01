
use std::str::from_utf8;

use needletail::parser::SequenceRecord;

use crate::utils::CompressionExt;
use crate::tools::error::ToolError;

pub fn get_compression_writer(
    output: &std::path::PathBuf,
    output_format: &Option<niffler::compression::Format>, 
    compression_level: &Option<niffler::compression::Level>
) -> Result<Box<dyn std::io::Write>, ToolError> {

    let file_handle = std::io::BufWriter::new(std::fs::File::create(&output)?);

    let fmt = match output_format {
        None => niffler::Format::from_path(&output),
        Some(format) => *format,
    };

    let level = match compression_level {
        Some(level) => *level,
        None => niffler::compression::Level::Six
    };
    
    Ok(niffler::get_writer(Box::new(file_handle), fmt, level)?)

}

pub fn get_seq_record_identifier(rec: &SequenceRecord) -> Result<String, ToolError> {
    match from_utf8(rec.id())?.split(' ').next() {
        Some(rec_id) => Ok(rec_id.to_string()),
        None => return Err(ToolError::RecordIdentifierNotParsed)
    }
}

pub fn get_seq_record_identifier_description(rec: &SequenceRecord) -> Result<(String, String), ToolError> {
    let mut header = from_utf8(rec.id())?.split(' ');
    let rec_id = match header.next() {
        Some(rec_id) => rec_id.to_string(),
        None => return Err(ToolError::RecordIdentifierNotParsed)
    };
    let rec_descr = match header.next() {
        Some(rec_descr) => rec_descr.to_string(),
        None => String::from("")
    };
    Ok((rec_id, rec_descr))
}