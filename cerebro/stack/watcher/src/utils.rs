use std::collections::HashMap;
use std::{io::Write, time::Duration};
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use cerebro_fs::client::UploadConfig;
use cerebro_model::api::files::model::{WatcherFormat, WatcherConfig};
use env_logger::Builder;
use env_logger::fmt::Color;
use log::{LevelFilter, Level};
use uuid::Uuid;

use crate::error::WatcherError;
use crate::terminal::WatchArgs;

pub const _CRATE_VERSION: &'static str = env!("CARGO_PKG_VERSION", "Failed to get the crate version at compile time - this is not good!");

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
            Some(Some("lzma")) => Self::Lzma,
            _ => Self::No,
        }
    }
}

pub trait StringUtils {
    fn substring(&self, start: usize, len: usize) -> Self;
}

impl StringUtils for String {
    fn substring(&self, start: usize, len: usize) -> Self {
        self.chars().skip(start).take(len).collect()
    }
}


pub trait UuidUtils {
    fn shorten(&self, len: usize) -> String;
}

impl UuidUtils for uuid::Uuid {
    fn shorten(&self, len: usize) -> String {
        self.to_string().substring(0, len)
    }
}

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

pub trait WatcherConfigArgs {
    fn from_args(watch_args: &WatchArgs) -> WatcherConfig;
}

impl WatcherConfigArgs for WatcherConfig {
    fn from_args(watch_args: &WatchArgs) -> WatcherConfig {
        WatcherConfig { 
            id: Uuid::new_v4().to_string(), 
            name: watch_args.name.clone(), 
            location: watch_args.location.clone(), 
            team_name: watch_args.team_name.clone(), 
            db_name: watch_args.db_name.clone(), 
            format: watch_args.format.clone(),
            interval: Duration::from_secs(watch_args.interval), 
            timeout: Duration::from_secs(watch_args.timeout), 
            timeout_interval: Duration::from_secs(watch_args.timeout_interval)
        }
    }
}

pub trait UploadConfigArgs {
    fn from_args(watch_args: &WatchArgs) -> UploadConfig;
}


impl UploadConfigArgs for UploadConfig {
    fn from_args(watch_args: &WatchArgs) -> UploadConfig {
        Self {
            data_center: watch_args.data_center.clone(),
            replication: watch_args.replication.clone(),
            ..Default::default()
        }
    }
}

pub trait CerebroClientConfigArgs {
    fn from_args(app_args: &crate::terminal::App) -> crate::watcher::CerebroClientConfig;
}

impl CerebroClientConfigArgs for crate::watcher::CerebroClientConfig {
    fn from_args(app_args: &crate::terminal::App) -> crate::watcher::CerebroClientConfig {
        crate::watcher::CerebroClientConfig { 
            api_url: app_args.url.clone(),
            api_token: app_args.token.clone(),
            api_token_file: app_args.token_file.clone(),
            _danger_invalid_certificate: app_args.danger_invalid_certificate,
            fs_url: app_args.fs_url.clone(),
            fs_port: app_args.fs_port.clone()
        }
    }
}

pub trait FileGetter {
    fn get_fastq_dir(&self, input_path: &PathBuf) -> Result<PathBuf, WatcherError>;
    fn get_fastq_files(&self, dir: &PathBuf, glob: Option<String>) -> Result<HashMap<String, Vec<PathBuf>>, WatcherError>;
}

impl FileGetter for WatcherFormat {
    fn get_fastq_dir(&self, input_path: &PathBuf) -> Result<PathBuf, WatcherError> {
        match self {
            WatcherFormat::Fastq => Ok(input_path.to_path_buf()),
            WatcherFormat::Iseq => {
                let alignment_path = input_path.join("Alignment_1");
                let latest = find_latest_alignment_directory(&alignment_path)?;
                let fastq_path = alignment_path.join(latest).join("Fastq");

                if fastq_path.exists() && fastq_path.is_dir() {
                    Ok(fastq_path)
                } else {
                    Err(WatcherError::InvalidDirectoryStructure(fastq_path))
                }
            }
            WatcherFormat::Nextseq => {
                let analysis_path = input_path.join("Analysis");
                let latest = find_latest_analysis_directory(&analysis_path)?;
                let fastq_path = analysis_path.join(latest).join("Data").join("fastq");

                if fastq_path.exists() && fastq_path.is_dir() {
                    Ok(fastq_path)
                } else {
                    Err(WatcherError::InvalidDirectoryStructure(fastq_path))
                }
            }
        }
    }
    fn get_fastq_files(&self, input_path: &PathBuf, glob: Option<String>) -> Result<HashMap<String, Vec<PathBuf>>, WatcherError> {
        
        let fastq_dir = self.get_fastq_dir(input_path)?;

        match self {
            WatcherFormat::Fastq => {
                get_read_files(&fastq_dir, &glob.unwrap_or("*.fastq.gz".to_string()), false)
            },
            WatcherFormat::Iseq => {
                get_read_files(&fastq_dir, "*_{L001_R1_001,L001_R2_001}.fastq.gz", false)
            },
            WatcherFormat::Nextseq => {
                get_read_files(&fastq_dir, "*_{R1_001,R2_001}.fastq.gz", false)
            }
        }

    }
}

fn find_latest_analysis_directory(path: &Path) -> Result<PathBuf, WatcherError> {
    let mut latest_number = None;
    let mut latest_dir = None;

    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            if let Some(dir_name) = entry.file_name().to_str() {
                if let Ok(number) = dir_name.parse::<u32>() {
                    if latest_number.is_none() || number > latest_number.unwrap() {
                        latest_number = Some(number);
                        latest_dir = Some(entry.path());
                    }
                }
            }
        }
    }

    latest_dir.ok_or_else(|| WatcherError::InvalidLatestAnalysisDirectory(path.to_path_buf()))
}

fn find_latest_alignment_directory(path: &Path) -> Result<PathBuf, WatcherError> {
    let mut latest_dir_name: Option<String> = None;
    let mut latest_dir: Option<PathBuf> = None;

    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            if let Some(dir_name) = entry.file_name().to_str() {
                if latest_dir_name.is_none() || dir_name > latest_dir_name.as_deref().unwrap() {
                    latest_dir_name = Some(dir_name.to_string());
                    latest_dir = Some(entry.path());
                }
            }
        }
    }

    latest_dir.ok_or_else(|| WatcherError::InvalidLatestAlignmentDirectory(path.to_path_buf()))
}



fn to_lexical_absolute(path: &PathBuf) -> Result<PathBuf, WatcherError> {
    let mut absolute = if path.is_absolute() {
        PathBuf::new()
    } else {
        std::env::current_dir()?
    };
    for component in path.components() {
        match component {
            std::path::Component::CurDir => {},
            std::path::Component::ParentDir => { absolute.pop(); },
            component @ _ => absolute.push(component.as_os_str()),
        }
    }
    Ok(absolute)
}


// A helper function to get read files from a suitable glob match
pub fn get_read_files(directory: &Path, glob: &str, symlinks: bool) -> Result<HashMap<String, Vec<PathBuf>>, WatcherError> {
   
    let glob = wax::Glob::new(glob).map_err(|_| WatcherError::GlobCreate(glob.to_string()))?;

    // Get potentially paired file paths into a HashMap
    let mut files = HashMap::new();
    for entry in glob.walk_with_behavior(directory, match symlinks { true => wax::LinkBehavior::ReadTarget, false => wax::LinkBehavior::ReadFile }) {
        let entry = entry.map_err(|_| WatcherError::GlobWalk(directory.to_path_buf()))?;
        let file_path = match symlinks {
            true => entry.path().canonicalize()?,
            false => to_lexical_absolute(&entry.path().to_path_buf())?
        };
        let sample_id = entry.matched().get(1).ok_or_else(|| WatcherError::GlobMatchSampleIdentifier(file_path.to_path_buf()))?;
        log::debug!("Sample sheet utility - [{:?}] - detected paired-end file: {:?}", sample_id, file_path);
        files.entry(sample_id.to_owned()).or_insert_with(Vec::new).push(file_path.to_path_buf());
    } 

    // For paired files the entries are unsorted! Sort them here by their names 
    // this will only work for traditional reverse/forward names like R1 and R2

    let mut sorted = HashMap::new();
    for (sample_id, mut file_paths) in files.clone().into_iter() {
        file_paths.sort();
        sorted.insert(sample_id, file_paths);
    }
    
    log::info!("{:#?}", files);

    Ok(sorted)
}