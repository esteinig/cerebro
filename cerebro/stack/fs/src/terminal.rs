use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

use cerebro_model::api::files::model::FileType;
use cerebro_model::api::files::retention::{RetentionClass, StorageTier};
use cerebro_client::terminal::StatusArgs;
use cerebro_client::terminal::LoginArgs;

/// Cerebro: file system anmd storage operations
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-fs", version)]
pub struct App {
    /// API URL
    #[clap(
        long, 
        short = 'u', 
        default_value = "http://localhost:8080", 
        env = "CEREBRO_API_URL"
    )]
    pub url: String,
    /// API token - usually provided with CEREBRO_API_TOKEN
    #[clap(
        long, 
        short = 'e', 
        env = "CEREBRO_API_TOKEN",
        hide_env_values = true
    )]
    pub token: Option<String>,
    /// API token file - can be set from environment variable
    #[clap(
        long, 
        short = 'f', 
        env = "CEREBRO_API_TOKEN_FILE"
    )]
    pub token_file: Option<PathBuf>,
    /// User team name or identifier for requests that require team specification 
    #[clap(
        long, 
        short = 't', 
        env = "CEREBRO_USER_TEAM",
        hide_env_values = true
    )]
    pub team: Option<String>,
    /// Team database name or identifier for requests that require database access 
    #[clap(
        long, 
        short = 'd', 
        env = "CEREBRO_USER_DB",
        hide_env_values = true
    )]
    pub db: Option<String>,
    /// Team database project name or identifier for requests that require project access 
    #[clap(
        long, 
        short = 'p', 
        env = "CEREBRO_USER_PROJECT",
        hide_env_values = true
    )]
    pub project: Option<String>,
    /// SeaweedFS master node address
    #[clap(
        long, 
        short = 'a', 
        default_value = "http://localhost", 
        env = "CEREBRO_FS_URL"
    )]
    pub fs_url: String,
    /// SeaweedFS master node port
    #[clap(
        long, 
        short = 'm',
        env = "CEREBRO_FS_PORT",
        default_value = "9333", 
    )]
    pub fs_port: String,
    /// SeaweedFS filer HTTP API base URL (path-addressed access)
    #[clap(
        long,
        env = "CEREBRO_FS_FILER_URL",
        default_value = "http://localhost:8888"
    )]
    pub fs_filer_url: String,
    /// Object I/O access mode: weed (fid-addressed, default) or filer (path-addressed)
    #[clap(
        long,
        value_enum,
        env = "CEREBRO_FS_ACCESS",
        default_value = "weed"
    )]
    pub fs_access: crate::config::FsAccessMode,
    /// SSL certificate verification is ignored [DANGER]
    #[clap(
        long, 
        env = "CEREBRO_DANGER_ACCEPT_INVALID_TLS_CERTIFICATE"
    )]
    pub danger_invalid_certificate: bool,

    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    // Status of Cerebro API
    Ping(StatusArgs),
    // Login to Cerebro API
    Login(LoginArgs),
    /// Upload files to CerebroFS and register files with CerebroAPI
    Upload(UploadFileArgs),
    /// Download of files from CerebroFS
    Download(DownloadFileArgs),
    /// Report and initiate archival (Glacier) restores for cold-tier files
    Restore(RestoreFileArgs),
    /// Verify file integrity against registered BLAKE3 hashes (optionally repair)
    Verify(VerifyFileArgs),
    /// Report health of the Cerebro FS topology (master, filer)
    Health,
    /// Delete a file from CerebroFS
    Delete(DeleteFileArgs),
    /// List accessible files from CerebroFS
    List(ListFileArgs),
    /// Stage samples in CerebroAPI / CerebroFS for production pipelines
    Stage(StageFileArgs),
    /// Get the SeaweedFS executable
    Weed(GetWeedArgs),
}

#[derive(Debug, Args)]
pub struct StageFileArgs {
    /// Staged sample model (.json)
    #[clap(long, short = 'j')]
    pub json: PathBuf,
    /// Stage file directory
    #[clap(long, short = 'o', default_value=".")]
    pub outdir: PathBuf,
    /// Stage a file that contains the requested pipeline
    #[clap(long, short = 'p')]
    pub pipeline: Option<PathBuf>,
}


#[derive(Debug, Args)]
pub struct GetWeedArgs {
    /// Executable output directory
    #[clap(long, short = 'o', default_value=".")]
    pub outdir: PathBuf,
    /// Executable output directory
    #[clap(long, short = 'v', default_value="latest")]
    pub version: String,
}

#[derive(Debug, Args)]
pub struct DeleteFileArgs {
    /// Files identifiers to delete (CerebroFS)
    #[clap(long, short = 'f', num_args(0..))]
    pub file_ids: Vec<String>,
    /// Sequence run identifier
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Sample identifier
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,
    /// Delete all files (requires confirmation)
    #[clap(long, short = 'a')]
    pub all: bool,
}


#[derive(Debug, Args)]
pub struct UploadFileArgs {
    /// Files to register
    #[clap(long, short = 'f', num_args(0..))]
    pub files: Vec<PathBuf>,
    /// File type
    #[clap(long, short = 't')]
    pub file_type: Option<FileType>,
    /// Sequence run identifier
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Biological sample identifier
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,
    /// Pipeline run identifier
    #[clap(long, short = 'p')]
    pub pipeline_id: Option<String>,
    /// File description
    #[clap(long, short = 'd')]
    pub description: Option<String>,
    /// Storage tier to register the file under
    #[clap(long, value_enum, default_value = "hot")]
    pub tier: StorageTier,
    /// Retention category to assign (resolved to an expiry by the retention policy)
    #[clap(long, value_enum, default_value = "diagnostic")]
    pub retention: RetentionClass,
    /// Register the file under legal hold (exempt from expiry)
    #[clap(long)]
    pub legal_hold: bool,
}



#[derive(Debug, Args)]
pub struct DownloadFileArgs {
    /// File identifiers to download directly (SeaweedFS fid or filer path).
    /// Bypasses the API lookup; integrity verification is unavailable for these.
    #[clap(long, short = 'f', num_args(0..))]
    pub fids: Vec<String>,
    /// Sequence run identifier (lists and downloads all registered files for the run)
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Restrict a run download to a single biological sample identifier
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,
    /// Output directory for downloaded files
    #[clap(long, short = 'o', default_value = ".")]
    pub outdir: PathBuf,
    /// Verify each downloaded file against its registered BLAKE3 hash
    #[clap(long)]
    pub verify: bool,
}

#[derive(Debug, Args)]
pub struct RestoreFileArgs {
    /// Sequence run identifier to evaluate for archival restores
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Restrict to a single biological sample identifier
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,
}

#[derive(Debug, Args)]
pub struct VerifyFileArgs {
    /// Sequence run identifier to verify
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Restrict verification to a single biological sample identifier
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,
    /// Attempt to repair a hash mismatch from an alternate replica
    #[clap(long)]
    pub repair: bool,
}

#[derive(Debug, Args)]
pub struct ListFileArgs {   
    /// Sequence run identifier
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Watcher identifier
    #[clap(long, short = 'w')]
    pub watcher_id: Option<String>,
    /// Return page of files
    #[clap(long, short = 'p', default_value="0")]
    pub page: u32,
    /// Files per page
    #[clap(long, short = 'l', default_value="1000")]
    pub limit: u32,
}

#[derive(Debug, Args)]
pub struct GlobalOptions {
}


pub fn get_styles() -> clap::builder::Styles {
	clap::builder::Styles::styled()
		.header(
			anstyle::Style::new()
				.bold()
				.underline()
				.fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
		)
		.literal(
			anstyle::Style::new()
				.bold()
				.fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
		)
}
