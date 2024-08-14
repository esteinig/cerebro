use std::collections::HashMap;
use std::io::Write;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use cerebro_client::client::CerebroClient;
use cerebro_fs::client::UploadConfig;
use cerebro_model::api::watchers::model::{WatcherFormat, ProductionWatcher};
use cerebro_model::api::watchers::schema::RegisterWatcherSchema;
use env_logger::Builder;
use env_logger::fmt::Color;
use log::{LevelFilter, Level};

use crate::error::WatcherError;
use crate::terminal::WatchArgs;

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
    fn from_args(watch_args: &WatchArgs, api_client: &CerebroClient) -> Result<ProductionWatcher, WatcherError>;
}

impl WatcherConfigArgs for ProductionWatcher {
    fn from_args(watch_args: &WatchArgs, api_client: &CerebroClient) -> Result<ProductionWatcher, WatcherError> {

        let registered_provided = watch_args.id.is_some() || watch_args.json.is_some();

        let new_provided = watch_args.name.is_some()
            && watch_args.location.is_some()
            && watch_args.format.is_some();

        match (registered_provided, new_provided) {
            (true, false) => {
                log::info!("Registered watcher arguments provided, getting production watcher...");

                let watcher_id = match (watch_args.json.clone(), watch_args.id.clone()) {
                    (Some(path), _) => RegisterWatcherSchema::from_json(&path)?.id,
                    (None, Some(id)) => id.to_owned(),
                    (None, None) => return Err(WatcherError::WatcherIdentifierArgNotFound)
                };

                let watchers = api_client.list_watchers(
                     Some(watcher_id), 
                     false
                )?;

                // With the optional identifier, the call returns a single item
                Ok(watchers[0].clone())

            },
            (false, true) => {
                log::info!("New watcher arguments provided, registering watcher...");

                let schema = RegisterWatcherSchema::new(
                    &watch_args.name.clone().unwrap(),
                    &watch_args.location.clone().unwrap(),
                    watch_args.format.clone().unwrap(),
                    watch_args.glob.clone()
                );

                api_client.register_watcher(
                    &schema, 
                    false
                )?;
                
                Ok(ProductionWatcher::from_schema(&schema))
                
            },
            _ => return Err(WatcherError::InvalidWatcherConfigArgs)
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
            ttl: watch_args.ttl.clone(),
            ..Default::default()
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
            WatcherFormat::Fastq | WatcherFormat::FastqPe => Ok(input_path.to_path_buf()),
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
        get_read_files(&self.get_fastq_dir(input_path)?, &glob.unwrap_or(self.default_glob()), false)
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