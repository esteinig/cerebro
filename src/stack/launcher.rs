use std::path::PathBuf;
use thiserror::Error;
use crate::pipeline::sheet::SampleSheet;
use crate::utils::UuidUtils;

use crate::tools::modules::slack::{SlackMessenger, SlackMessage, SlackMessageSectionBlock, TextObject};

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
    /// Indicates failure to parse sample sheet
    #[error("failed to read sample sheet")]
    SampleSheetNotRead,
    /// Indicates failure to find entries in sample sheet
    #[error("failed to detect entries in sample sheet")]
    SampleSheetEmpty,

    /// Indicates failure to detect sample sheet
    #[error("failed to detect sample sheet")]
    SampleSheetNotFound,

    /// Indicates failure to detect fastq sub-directory
    #[error("failed to detect fastq sub-directory")]
    FastqDirectoryNotFound,
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
        let sample_sheet_path = self.input_path.join(self.run_id.clone()).with_extension("csv");

        validation.sample_sheet_detected.checked = true;
        match sample_sheet_path.exists() {
            false => {
                log::warn!("[{run_id}] Failed to detect sample sheet: {}", sample_sheet_path.display());

                // Send slack message and return error as downstream checks depend on sample sheet
                self.slack.send(&self.slack_message.input_validation(
                    &validation, &run_id, &self.launch_id,
                )).map_err(|_| WorkflowLauncherError::SlackMessageNotSent)?;

                return Err(WorkflowLauncherError::SampleSheetNotFound)  
            },
            true => {
                log::info!("[{run_id}] Sample sheet detected");
                validation.sample_sheet_detected.pass = true;
            }
        }

        // Does the read input folder exist?
        let fastq_subdir_path = self.input_path.join("fastq");

        validation.fastq_subdir_detected.checked = true;
        match fastq_subdir_path.exists() {
            false => {
                log::warn!("[{run_id}] Failed to detect fastq sub-directory: {}", fastq_subdir_path.display());

                // Send slack message and return error as downstream checks depend on sample sheet
                self.slack.send(&self.slack_message.input_validation(
                    &validation, &run_id, &self.launch_id,
                )).map_err(|_| WorkflowLauncherError::SlackMessageNotSent)?;

                return Err(WorkflowLauncherError::FastqDirectoryNotFound)  
            },
            true => {
                log::info!("[{run_id}] Fastq sub-directory detected");
                validation.fastq_subdir_detected.pass = true;
            }
        }

        // Sample Sheet checks:

        // Can the sample sheet be parsed? Includes check that sample identifiers are unique for this sample sheet. 

        // Failure may be due to malformatted CSV or missing required fields. We should output the error message into 
        // the slack message for this (since it can be opaque why initial parsing failed otherwise)

        validation.sample_sheet_parsed.checked = true;
        let sample_sheet = match SampleSheet::from(&sample_sheet_path) {
            Ok(sample_sheet) => {
                log::info!("[{run_id}] Sample sheet parsed successfully");
                validation.sample_sheet_parsed.pass = true;
                sample_sheet
            },
            Err(err) => {
                log::warn!("[{run_id}] Failed to parse sample sheet: {}", sample_sheet_path.display());

                // Full JSON style error format for additional information
                validation.sample_sheet_parsed.error_message = format!("```{}```", format!("{err:#?}").replace("`", ""));

                // Send slack message and return error as downstream checks depend on sample sheet
                self.slack.send(&self.slack_message.input_validation(
                    &validation, &run_id, &self.launch_id,
                )).map_err(|_| WorkflowLauncherError::SlackMessageNotSent)?;

                return Err(err).map_err(|_| WorkflowLauncherError::SampleSheetNotRead)  
            }
        };

        // Is the sample sheet empty?

        validation.sample_sheet_not_empty.checked = true;
        if !(sample_sheet.entries.len() > 0) {
            log::warn!("[{run_id}] Sample sheet is empty: {}", sample_sheet_path.display());

            // Send slack message and return error as downstream checks depend on sample sheet
            self.slack.send(&self.slack_message.input_validation(
                &validation, &run_id, &self.launch_id,
            )).map_err(|_| WorkflowLauncherError::SlackMessageNotSent)?;

            return Err(WorkflowLauncherError::SampleSheetEmpty)
        } else {
            log::info!("[{run_id}] Sample sheet is not empty");
            validation.sample_sheet_not_empty.pass = true;
        }

        // // Are there missing required values for samples? - done automatically by CsvReader



        // // Are there non-allowed characters in the fields?

        // sample_sheet.validate_field_characters();

        // // Is the project name provided? [if not provided in launcher configuration]

        // sample_sheet.validate_project_name();

        // // Do the sample identifiers correspond to tagged format? [strict]

        // sample_sheet.validate_sample_name_tags();

        // // Is the run date format correct? [strict]

        // sample_sheet.validate_date_format();

        // Fastq checks:

        // Are the fastq files as specified by the sample sheet identifiers present?

        // Are the fastq files not empty?


        


        // Create validation summary and compose slack pass/fail message

        self.slack.send(&self.slack_message.input_validation(
            &validation, &run_id, &self.launch_id
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
        SlackMessage::new(&self.channel, &format!("New run input detected [{run_id}::{}]", launch_id.shorten(8)))
    }
    pub fn input_validation(&self, validation: &InputValidation, run_id: &str, launch_id: &uuid::Uuid) -> SlackMessage {
        match validation.pass() {
            true => SlackMessage::new(
                &self.channel, 
                &format!("Input validation passed [{run_id}::{}]", launch_id.shorten(8))
            ),
            false => {
                let msg = &format!("*Input validation failed* [{run_id}::{}]", launch_id.shorten(8));
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
        self.sample_sheet_detected.pass &&
        self.fastq_subdir_detected.pass &&
        self.sample_sheet_parsed.pass && 
        self.sample_sheet_not_empty.pass
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
