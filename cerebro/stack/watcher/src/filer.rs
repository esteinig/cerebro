use std::{collections::HashMap, path::{Path, PathBuf}};
use cerebro_client::client::CerebroClient;
use cerebro_fs::client::FileSystemClient;
use thiserror::Error;

use crate::watcher::SlackConfig;
use cerebro_model::api::files::model::WatcherConfig;
use crate::slack::{SlackMessage, SlackMessageSectionBlock, SlackMessenger, TextObject};
use cerebro_fs::client::UploadConfig;

#[derive(Error, Debug)]
pub enum WatchFilerError {
    /// Indicates failure input/output file
    #[error("failed to read file")]
    FileIO(#[from] std::io::Error),
    /// Indicates failure obtain path name
    #[error("failed to parse base name of: {0}")]
    PathBaseName(String),
    /// Indicates failure to send Slack message
    #[error("failed to send slack message")]
    SlackMessageNotSent,
    /// Indicates failure to validate inputs
    #[error("failed input validation")]
    InputValidationFailed,
    /// Indicates failure to parse sample sheet
    #[error("failed to read sample sheet")]
    _SampleSheetNotRead,
    /// Indicates failure to find entries in sample sheet
    #[error("failed to detect entries in sample sheet")]
    _SampleSheetEmpty,
    /// Indicates failure to detect sample sheet
    #[error("failed to detect sample sheet")]
    SampleSheetNotFound,
    /// Indicates failure to detect fastq sub-directory
    #[error("failed to detect fastq sub-directory")]
    FastqDirectoryNotFound,

    #[error("error in client lib")]
    HttpClientError(#[from] cerebro_client::error::HttpClientError),
    #[error("error in cerebro fs lib")]
    CerebroFsError(#[from] cerebro_fs::error::FileSystemError),

    /// Represents a failure to extract the sample identifier from a glob matched pattern
    #[error("failed to extract pattern match for sample identifier from: {0}")]
    GlobMatchSampleIdentifier(String),
    /// Represents a failure to 
    #[error("failed to extract an entry from the globbed walk through directory: {0}")]
    GlobWalk(String),
    /// Represents a failure to 
    #[error("failed to create a glob matcher for pattern: {0}")]
    GlobCreate(String),
    /// Represents a failure to 
    #[error("failed to find paired files for sample: {0}")]
    GlobPairedFiles(String),
    /// Represents a failure to 
    #[error("failed to get inernal fastq path - has it been validated?")]
    FastqDirNotValidated,
}


pub struct FilerConfig {
    fastq_subdir: PathBuf,
    sample_sheet: Option<PathBuf>
}

impl Default for FilerConfig {
    fn default() -> Self {
        Self {
            fastq_subdir: PathBuf::from("fastq"),
            sample_sheet: None
        }
    }
}


pub struct SlackTools {
    client: SlackMessenger,
    message: SlackMessageGenerator
}

#[derive(Clone, Debug)]
pub struct CerebroClientConfig {
    api_url: String,
    api_token: Option<String>,
    api_token_file: Option<PathBuf>,
    _danger_invalid_certificate: bool,
    fs_url: String,
    fs_port: String,
}
impl Default for CerebroClientConfig {
    fn default() -> Self {
        Self {
            api_url: String::from("http://api.dev.cerebro.localhost"),
            api_token: std::env::var("CEREBRO_API_TOKEN").ok(),
            api_token_file: None,
            _danger_invalid_certificate: false,
            fs_url: String::from("http://fs.dev.cerebro.localhost"),
            fs_port: String::from("9333"),
        }
    }
}

/// Input validation and workflow launcher for production
pub struct WatchFiler {
    pub base_path: PathBuf,
    pub input_path: PathBuf,
    pub slack: Option<SlackTools>,
    pub run_id: String,
    pub fastq_path: Option<PathBuf>,
    pub watcher_config: WatcherConfig,
    pub filer_config: FilerConfig,
    pub fs_client: FileSystemClient,
}

impl WatchFiler {
    pub fn new(base_path: PathBuf, input_path: PathBuf, watcher_config: WatcherConfig, cerebro_client_config: CerebroClientConfig, slack_config: Option<SlackConfig>,) -> Result<Self, WatchFilerError> {
        
        // Input path directory is the run identifier
        let path_string = input_path.display().to_string();

        let run_id = match input_path.file_name() {
            Some(name) => match name.to_os_string().into_string() {
                Ok(run_id) => run_id, 
                Err(_) => return Err(WatchFilerError::PathBaseName(path_string))
            },
            None => return Err(WatchFilerError::PathBaseName(path_string))
        };

       let slack = match slack_config {
            Some(slack_config) => {
                // Slack notification setup
                let slack = SlackMessenger::new(&slack_config.token);
                let slack_message = SlackMessageGenerator::new(slack_config.channel);

                // Send initialisation message
                slack.send(
                    &slack_message.input_detected(&run_id, &watcher_config.name, &watcher_config.location)
                ).map_err(|_| WatchFilerError::SlackMessageNotSent)?;

                Some(SlackTools {
                    client: slack,
                    message: slack_message
                })

            },
            _ => None
        };

        let api_client = CerebroClient::new(
            &cerebro_client_config.api_url,
            &cerebro_client_config.api_token,
            false,
            false,
            &cerebro_client_config.api_token_file
        )?;

        let fs_client = FileSystemClient::new(
            &api_client, 
            &cerebro_client_config.fs_url, 
            &cerebro_client_config.fs_port
        );

        log::info!("Checking status of Cerebro API at {}",  &cerebro_client_config.api_url);
        api_client.ping_servers()?;

        log::info!("Checking status of SeaweedFS master at {}",  &cerebro_client_config.fs_url);
        fs_client.ping_status()?;

        
        Ok(Self { base_path, input_path, run_id, slack, watcher_config, filer_config: FilerConfig::default(), fs_client, fastq_path: None })
    }

    /// Upload and register target files to CerebroFS
    pub fn upload_and_register_illumina_pe(
        &self,
        team_name: &str,
        db_name: &str,
        paired_glob: &str, 
        allow_single: bool, 
        follow_symlinks: bool, 
    )-> Result<(), WatchFilerError> {

        let run_id = self.run_id.clone();
        let watcher_name = self.watcher_config.name.clone();
        let watcher_loc = self.watcher_config.location.clone();

        let dir = match &self.fastq_path {
            Some(path) => path, None => return Err(WatchFilerError::FastqDirNotValidated)
        };

        let files = get_paired_files(
            &dir, paired_glob, allow_single, follow_symlinks
        )?;

        if files.is_empty() {
            log::warn!("[{watcher_name}@{watcher_loc}::{run_id}] FAILED to detect Illumina PE read files with extension glob `{paired_glob}` in {}", dir.display());
            return Ok(())
        }

        log::info!("[{watcher_name}@{watcher_loc}::{run_id}] Uploading and registering files with Cerebro FS API");

        log::info!("Starting file processing and upload...");
        self.fs_client.upload_illumina_pe(
            &files,  
            team_name,
            db_name,
            UploadConfig::default(),
            WatcherConfig::default()
        )?;
        
        Ok(())
    }

    /// Validates presence and integrity of input sample sheet and read files 
    pub fn validate_inputs(&mut self) -> Result<(), WatchFilerError> {

        let run_id = self.run_id.clone();
        let watcher_name = self.watcher_config.name.clone();
        let watcher_loc = self.watcher_config.location.clone();

        // Conduct all checks first
        let mut validation = InputValidation::new();

        match &self.filer_config.sample_sheet {
            Some(sample_sheet) => {
                match sample_sheet.exists() {
                    false => {
                        log::warn!("[{watcher_name}@{watcher_loc}::{run_id}] Failed to detect sample sheet: {}", sample_sheet.display());
                        validation.sample_sheet_detected.checked = true;
                        validation.sample_sheet_detected.pass = false;

                    },
                    true => {
                        log::info!("[{watcher_name}@{watcher_loc}::{run_id} Sample sheet detected");
                        validation.sample_sheet_detected.checked = true;
                        validation.sample_sheet_detected.pass = true;
                    }
                }
            },
            None => {
                // Sample sheet is optional
                validation.sample_sheet_detected.checked = false;
                validation.sample_sheet_detected.pass = false;
            }
        };

        // Does the read input folder exist?
        let fastq_subdir_path = self.input_path.join(
            self.filer_config.fastq_subdir.clone()
        );

        match fastq_subdir_path.exists() {
            false => {
                log::warn!("[{watcher_name}@{watcher_loc}::{run_id}] Failed to detect fastq sub-directory: {}", fastq_subdir_path.display());

                validation.fastq_subdir_detected.checked = true;
                validation.fastq_subdir_detected.pass = false;
            },
            true => {
                log::info!("[{watcher_name}@{watcher_loc}::{run_id}] Fastq sub-directory detected");

                validation.fastq_subdir_detected.checked = true;
                validation.fastq_subdir_detected.pass = true;
            }
        }


        // Create validation summary and compose slack pass/fail message
        if let Some(slack_tools) = &self.slack {
            slack_tools.client.send(&slack_tools.message.input_validation(
                &validation, &run_id, &self.watcher_config.name, &self.watcher_config.location
            )).map_err(|_| WatchFilerError::SlackMessageNotSent)?;
        }

        // Validate on InputValidation scheme
        if validation.pass() {
            self.fastq_path = Some(fastq_subdir_path);
            Ok(())
        } else {
            Err(WatchFilerError::InputValidationFailed)
        }
    }
}

pub struct SlackMessageGenerator {
    pub channel: String
}
impl SlackMessageGenerator {
    pub fn new(channel: String) -> Self {
        Self { channel }
    }
    pub fn input_detected(&self, run_id: &str, watcher_name: &str, watcher_loc: &str) -> SlackMessage {
        SlackMessage::new(&self.channel, &format!("[{watcher_name}@{watcher_loc}::{run_id} New run input detected"))
    }
    pub fn input_validation(&self, validation: &InputValidation, run_id: &str, watcher_name: &str, watcher_loc: &str) -> SlackMessage {
        match validation.pass() {
            true => SlackMessage::new(
                &self.channel, 
                &format!("[{watcher_name}@{watcher_loc}::{run_id}] Input validation passed")
            ),
            false => {
                let msg = &format!("[{watcher_name}@{watcher_loc}::{run_id}] *Input validation failed*");
                SlackMessage::from(&self.channel, vec![
                    vec![SlackMessageSectionBlock::new(TextObject::markdown(msg))],
                    validation.get_message_blocks(&run_id)
                ].concat())
            }
        }
    }
}

pub struct InputValidationChecks {
    pub pass: bool,
    pub checked: bool,
    pub error_message: String
}
impl Default for InputValidationChecks {
    fn default() -> Self {
        Self { pass: false, checked: false, error_message: String::from("") }
    }
}

pub struct InputValidation {
    pub sample_sheet_detected: InputValidationChecks,
    pub fastq_subdir_detected: InputValidationChecks,
    pub sample_sheet_parsed: InputValidationChecks,
    pub sample_sheet_not_empty: InputValidationChecks,
}
impl InputValidation {
    pub fn new() -> Self {
        Self { 
            sample_sheet_detected: Default::default(), 
            fastq_subdir_detected:  Default::default(),
            sample_sheet_parsed:  Default::default(),
            sample_sheet_not_empty:  Default::default()
        }
    }
    pub fn pass(&self) -> bool {
        self.fastq_subdir_detected.pass
    }
    pub fn get_message_blocks(&self, run_id: &str) -> Vec<SlackMessageSectionBlock> {
        vec![
            self.get_block(&self.sample_sheet_detected, &format!("Is the sample sheet `{run_id}.csv` present?")),
            self.get_block(&self.fastq_subdir_detected, "Is the sub-directory `fastq` present?"),
            self.get_block(&self.sample_sheet_parsed, "Is the sample sheet formatted correctly?"),
            self.get_block(&self.sample_sheet_not_empty, "Is the sample sheet not empty?")
        ]
    }
    pub fn get_block(&self, validation: &InputValidationChecks, msg: &str)-> SlackMessageSectionBlock {
        SlackMessageSectionBlock::new(TextObject::markdown(
        &format!("{} >> {msg} \n {}", 
                match validation.pass { true => "✅", false => match validation.checked { true => "❌", false => "✖" }},
                match validation.error_message.is_empty() { true => String::new(), false => format!("{}", validation.error_message)}
            )
        ))
    }
}


// A helper function to get paired files from a suitable glob match
pub fn get_paired_files(directory: &Path, paired_glob: &str, single: bool, symlinks: bool) -> Result<HashMap<String, Vec<PathBuf>>, WatchFilerError> {
   
    let glob = wax::Glob::new(paired_glob).map_err(|_| WatchFilerError::GlobCreate(paired_glob.to_string()))?;

    // Get potentially paired file paths into a HashMap
    let mut paired_files = HashMap::new();
    for entry in glob.walk_with_behavior(directory, match symlinks { true => wax::LinkBehavior::ReadTarget, false => wax::LinkBehavior::ReadFile }) {
        let entry = entry.map_err(|_| WatchFilerError::GlobWalk(format!("{:?}", &directory)))?;
        let file_path = match symlinks {
            true => entry.path().canonicalize()?,
            false => to_lexical_absolute(&entry.path().to_path_buf())?
        };
        let sample_id = entry.matched().get(1).ok_or_else(||WatchFilerError::GlobMatchSampleIdentifier(format!("{:?}", file_path)))?;
        log::debug!("Sample sheet utility - [{:?}] - detected paired-end file: {:?}", sample_id, file_path);
        paired_files.entry(sample_id.to_owned()).or_insert_with(Vec::new).push(file_path.to_path_buf());
    } 

    // Check the entries of the HashMap if single files are not allowed:
    let paired_files = match single {
        false => {
            let mut sorted = HashMap::new();
            for (sample_id, mut file_paths) in paired_files.clone().into_iter() {
                if file_paths.len() != 2 {
                    return Err(WatchFilerError::GlobPairedFiles(sample_id.to_string()))
                }
                // For paired files the entries are unsorted! Sort them here by their names 
                // this will only work for traditional reverse/forward names like R1 and R2
                file_paths.sort();
                sorted.insert(sample_id, file_paths);
            }
            sorted
        },
        true => paired_files
    };
    
    log::info!("{:#?}", paired_files);

    Ok(paired_files)
}


fn to_lexical_absolute(path: &PathBuf) -> std::io::Result<PathBuf> {
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