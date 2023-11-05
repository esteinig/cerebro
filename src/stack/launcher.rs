use std::path::PathBuf;
use thiserror::Error;
use crate::utils::UuidUtils;

use crate::tools::modules::slack::{SlackMessenger, SlackMessage};

#[derive(Error, Debug)]
pub enum WorkflowLauncherError {
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
}

/// Input validation and workflow launcher for production
pub struct WorkflowLauncher {
    pub base_path: PathBuf,
    pub input_path: PathBuf,
    pub slack: SlackMessenger,
    pub slack_message: SlackMessageGenerator,
    pub run_id: String,
    pub launch_id: uuid::Uuid
}

impl WorkflowLauncher {
    pub fn new(base_path: PathBuf, input_path: PathBuf, slack_channel: String, slack_token: String) -> Result<Self, WorkflowLauncherError> {
        
        // Input path directory is the run identifier
        let path_string = input_path.display().to_string();

        let run_id = match input_path.file_name() {
            Some(name) => match name.to_os_string().into_string() {
                Ok(run_id) => run_id, 
                Err(_) => return Err(WorkflowLauncherError::PathBaseName(path_string))
            },
            None => return Err(WorkflowLauncherError::PathBaseName(path_string))
        };

        // Launch identifier - used in the execution directory
        let launch_id = uuid::Uuid::new_v4();

        // Slack notification setup
        let slack = SlackMessenger::new(&slack_token);
        let slack_message = SlackMessageGenerator::new(slack_channel);

        // Send initialisation message
        slack.send(
            &slack_message.input_detected(&run_id, &launch_id)
        ).map_err(|_| WorkflowLauncherError::SlackMessageNotSent)?;

        Ok(Self { base_path, input_path, run_id, slack_message, slack, launch_id })
    }
    /// Validates presence and integrity of input sample sheet and read files 
    pub fn validate_inputs(&self) -> Result<(), WorkflowLauncherError> {

        let run_id = self.run_id.clone();

        // Conduct all checks first
        let mut validation = InputValidation::new();

        // Does the sample sheet {run_id}.csv exist?
        let sample_sheet = self.input_path.join(self.run_id.clone()).with_extension("csv");

        match sample_sheet.exists() {
            false => {
                log::warn!("[{run_id}] Failed to detect sample sheet: {}", sample_sheet.display());
            },
            true => {
                log::info!("[{run_id}] Sample sheet detected ✅");
                validation.sample_sheet_detected = true;
            }
        }

        // Does the read input folder exist?
        let fastq_subdir = self.input_path.join("fastq");

        match fastq_subdir.exists() {
            false => {
                log::warn!("[{run_id}] Failed to detect fastq sub-directory: {}", fastq_subdir.display());
            },
            true => {
                log::info!("[{run_id}] Fastq sub-directory detected ✅");
                validation.fastq_subdir_detected = true;
            }
        }

        // Fastq checks:

        // Are the fastq files as specified by the sample sheet identifiers present?

        // Are the fastq files not empty?


        
        // Sample Sheet checks:

        // Can the sample sheet be parsed?

        // Are there missing required values for samples?

        // Are there non-allowed characters in the fields?

        // Is the project name provided? [if not provided in launcher configuration]

        // Do the sample identifiers correspond to tagged format? [strict]

        // Is the run date format correct? [strict]


        // Create validation summary and compose slack pass/fail message

        self.slack.send(&self.slack_message.input_validation(
            validation.pass(), &run_id, &self.launch_id
        )).map_err(|_| WorkflowLauncherError::SlackMessageNotSent)?;

        if validation.pass() {
            Ok(())
        } else {
            Err(WorkflowLauncherError::InputValidationFailed)
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
    pub fn input_detected(&self, run_id: &str, launch_id: &uuid::Uuid) -> SlackMessage {
        SlackMessage::new(&self.channel, &format!("👀 New run input detected [run_id=*{run_id}* launch_id={}]", launch_id.shorten(8)))
    }
    pub fn input_validation(&self, pass: bool, run_id: &str, launch_id: &uuid::Uuid) -> SlackMessage {
        let msg = match pass {
            true => format!("✅ Input validation passed [run_id=*{run_id}* launch_id={}]", launch_id.shorten(8)),
            false => format!("❌ Input validation failed [run_id=*{run_id}* launch_id={}]", launch_id.shorten(8))
        };
        SlackMessage::new(&self.channel, &msg)
    }
}

pub struct InputValidation {
    pub sample_sheet_detected: bool,
    pub fastq_subdir_detected: bool,
}
impl InputValidation {
    pub fn new() -> Self {
        Self { 
            sample_sheet_detected: false, 
            fastq_subdir_detected: false 
        }
    }
    pub fn pass(&self) -> bool {
        self.sample_sheet_detected &&
        self.fastq_subdir_detected
    }
}
