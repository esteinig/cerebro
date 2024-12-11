use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

use cerebro_model::api::files::model::FileType;
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
        default_value = "http://api.cerebro.localhost", 
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
        default_value = "http://fs.cerebro.localhost", 
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
}



#[derive(Debug, Args)]
pub struct DownloadFileArgs {   

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
