use csv::{Reader, ReaderBuilder, Writer, WriterBuilder};
use env_logger::{fmt::Color, Builder};
use log::{Level, LevelFilter};
use needletail::{parse_fastx_file, FastxReader};
use niffler::get_writer;
use niffler::seek::compression::ReadSeek;
use niffler::seek::get_reader;
use serde::de::IntoDeserializer;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::io::{BufReader, BufWriter, Cursor, Read, Write};
use std::io::BufRead;
use crate::error::WorkflowError;

pub fn init_logger() {

    Builder::new()
        .format(|buf, record| {
            let timestamp = buf.timestamp();

            let mut red_style = buf.style();
            red_style.set_color(Color::Red).set_bold(true);
            let mut green_style = buf.style();
            green_style.set_color(Color::Green).set_bold(true);
            let mut white_style = buf.style();
            white_style.set_color(Color::White).set_bold(false);
            let mut orange_style = buf.style();
            orange_style.set_color(Color::Rgb(255, 102, 0)).set_bold(true);
            let mut apricot_style = buf.style();
            apricot_style.set_color(Color::Rgb(255, 195, 0)).set_bold(true);

            let msg = match record.level(){
                Level::Warn => (orange_style.value(record.level()), orange_style.value(record.args())),
                Level::Info => (green_style.value(record.level()), white_style.value(record.args())),
                Level::Debug => (apricot_style.value(record.level()), apricot_style.value(record.args())),
                Level::Error => (red_style.value(record.level()), red_style.value(record.args())),
                _ => (white_style.value(record.level()), white_style.value(record.args()))
            };

            writeln!(
                buf,
                "{} [{}] - {}",
                white_style.value(timestamp),
                msg.0,
                msg.1
            )
        })
        .filter(None, LevelFilter::Info)
        .init();
}


/// Enum to specify the type of file component to retrieve
pub enum FileComponent {
    /// The full file name including the extension
    FileName,
    /// The file name without the extension
    FileStem,
}

/// Extracts the specified file component from a `PathBuf` and returns it as a `String`.
///
/// # Arguments
///
/// * `path` - A `PathBuf` representing the file path.
/// * `component` - A `FileComponent` specifying whether to get the file name or the file stem.
///
/// # Returns
///
/// * `Result<String, DatabaseError>` - A `Result` containing the specified file component as a `String`
///   if successful, or a `DatabaseError` if an error occurs.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
/// use cipher::{get_file_component, FileComponent};
///
/// let path = PathBuf::from("/some/path/to/file.txt");
/// match get_file_component(path, FileComponent::FileName) {
///     Ok(file_name) => println!("File name: {}", file_name),
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
///
/// ```
/// use std::path::PathBuf;
/// use cipher::{get_file_component, FileComponent};
///
/// let path = PathBuf::from("/some/path/to/file.txt");
/// match get_file_component(path, FileComponent::FileStem) {
///     Ok(file_stem) => println!("File stem: {}", file_stem),
///     Err(e) => eprintln!("Error: {}", e),
/// }
/// ```
pub fn get_file_component(path: &PathBuf, component: FileComponent) -> Result<String, WorkflowError> {
    match component {
        FileComponent::FileName => {
            path.file_name()
                .ok_or(WorkflowError::FileNameConversionError)
                .and_then(|os_str| os_str.to_str().map(String::from).ok_or(WorkflowError::FileNameConversionError))
        }
        FileComponent::FileStem => {
            path.file_stem()
                .ok_or(WorkflowError::FileNameConversionError)
                .and_then(|os_str| os_str.to_str().map(String::from).ok_or(WorkflowError::FileNameConversionError))
        }
    }
}


pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

/// Attempts to infer the compression type from the file extension.
/// If the extension is not known, then Uncompressed is returned.
impl CompressionExt for niffler::compression::Format {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self {
        let path = Path::new(p);
        match path.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => Self::Gzip,
            Some(Some("bz") | Some("bz2")) => Self::Bzip,
            Some(Some("lzma") | Some("xz")) => Self::Lzma,
            _ => Self::No,
        }
    }
}

pub fn get_compression_writer(
    output: &std::path::PathBuf,
    output_format: &Option<niffler::compression::Format>, 
    compression_level: &Option<niffler::compression::Level>
) -> Result<Box<dyn std::io::Write>, WorkflowError> {

    let file_handle = std::io::BufWriter::new(std::fs::File::create(&output)?);

    let fmt = match output_format {
        None => niffler::Format::from_path(&output),
        Some(format) => *format,
    };

    let level = match compression_level {
        Some(level) => *level,
        None => niffler::compression::Level::Six
    };
    
    Ok(niffler::get_writer(
        Box::new(file_handle), 
        fmt, 
        level
    )?
    )

}


pub fn is_file_empty<P: AsRef<Path>>(path: P) -> Result<bool, WorkflowError> {
    let file = File::open(&path)?;
    
    // Use niffler to get a reader for the (possibly compressed) file
    let (mut reader, _format) = match get_reader(Box::new(file)) {
        Ok(reader_format) => reader_format,
        Err(niffler::Error::FileTooShort) => return Ok(true),
        Err(e) => return Err(WorkflowError::NifflerError(e)),
    };
    // Try to read the first byte
    let mut buffer = [0; 1];
    match reader.read(&mut buffer) {
        Ok(0) => Ok(true),
        Ok(_) => Ok(false), // Successfully read a byte, file is not empty
        Err(e) => Err(WorkflowError::IoError(e))
    }
}

pub fn parse_fastx_file_with_check<P: AsRef<Path>>(path: P) -> Result<Option<Box<dyn FastxReader>>, WorkflowError> {
    if is_file_empty(&path)? {
        Ok(None)
    } else {
        Ok(Some(parse_fastx_file(&path)?))
    }
}

pub fn get_file_by_name(path: &PathBuf, id: &str, extension: &str) -> Result<Option<PathBuf>, WorkflowError> {
    
    let file_path = path.join(format!("{id}{extension}"));
    
    match file_path.exists() && file_path.is_file() { 
        true => Ok(Some(file_path)), 
        false => Ok(None) 
    }
}
pub fn get_files_from_patterns(path: &PathBuf, patterns: &[&str]) -> Result<Option<Vec<PathBuf>>, WorkflowError> {
    
    // Need to have canonical path for GlobWalk
    let full_path = std::fs::canonicalize(path).map_err(|_| WorkflowError::InvalidReferencePath)?;

    let walker = globwalk::GlobWalkerBuilder::from_patterns(
        full_path, patterns
    )
        .max_depth(1)
        .follow_links(false)
        .build().map_err(|_| WorkflowError::GlobWalkBuilder)?
        .into_iter()
        .filter_map(Result::ok);
    
    let mut alignments = Vec::new();
    for file in walker {
        alignments.push(file.path().to_path_buf())
    };
    match alignments.is_empty() {
        true => Ok(None),
        false => Ok(Some(alignments))
    }
}


/// Utility function to extract the ID from a FASTQ record header.
///
/// # Arguments
///
/// * `id` - A byte slice containing the FASTQ record header.
///
/// # Returns
///
/// * `Result<String, WorkflowError>` - The extracted ID as a string on success, otherwise an error.
///
/// # Example
///
/// ```
/// let id = get_id(b"@read1 description").unwrap();
/// ```
pub fn get_id(id: &[u8]) -> Result<String, WorkflowError> {
    let header = std::str::from_utf8(id)?;
    let header_components = header
        .split_whitespace()
        .collect::<Vec<&str>>();
    
    if header_components.len() < 1 {
        return Err(WorkflowError::NeedletailFastqHeader)
    }
    let id = header_components[0].to_string();

    Ok(id)
}



pub fn get_tsv_reader(file: &Path, flexible: bool, header: bool) -> Result<Reader<Box<dyn ReadSeek>>, WorkflowError> {

    let buf_reader = BufReader::new(File::open(&file)?);
    let (reader, _format) = get_reader(Box::new(buf_reader))?;

    let tsv_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(header)
        .flexible(flexible) // Allows records with a different number of fields
        .from_reader(reader);

    Ok(tsv_reader)
}

pub fn get_tsv_writer(file: &Path, header: bool) -> Result<Writer<Box<dyn Write>>, WorkflowError> {
    
    let buf_writer = BufWriter::new(File::create(&file)?);
    let writer = get_writer(Box::new(buf_writer), niffler::Format::from_path(file), niffler::compression::Level::Six)?;

    let csv_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(header)
        .from_writer(writer);

    Ok(csv_writer)
}

pub fn write_tsv<T: Serialize>(data: &Vec<T>, file: &Path, header: bool) -> Result<(), WorkflowError> {

    let mut writer = get_tsv_writer(file, header)?;

    for value in data {
        // Serialize each value in the vector into the writer
        writer.serialize(&value)?;
    }

    // Flush and complete writing
    writer.flush()?;
    Ok(())
}

pub fn read_tsv<T: for<'de>Deserialize<'de>>(file: &Path, flexible: bool, header: bool) -> Result<Vec<T>, WorkflowError> {

    let mut reader = get_tsv_reader(file, flexible, header)?;

    let mut records = Vec::new();
    for record in reader.deserialize() {
        records.push(record?)
    }

    Ok(records)
}


pub fn read_tsv_skip<T: for<'de>Deserialize<'de>>(file: &Path, flexible: bool, header: bool, skip: char) -> Result<Vec<T>, WorkflowError> {

    let file = File::open(file)?;
    let reader = BufReader::new(file);

    // Filter lines starting with '#'
    let filtered_lines: Vec<String> = reader.lines()
        .filter_map(Result::ok)
        .filter(|line| !line.starts_with(skip))
        .collect();

    // Join filtered lines into a single String, separated by newlines
    let filtered_content = filtered_lines.join("\n");

    log::info!("Filtered content: {filtered_content}");

    // Use a Cursor to read this content as CSV input
    let mut csv_reader = ReaderBuilder::new()
        .has_headers(header)
        .flexible(flexible)
        .from_reader(Cursor::new(filtered_content));

    // Iterate over records
    let mut records = Vec::new();
    for result in csv_reader.deserialize() {
        let record = result?;
        records.push(record)
    }

    Ok(records)
}