use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

use cerebro_client::terminal::ApiStatusArgs;
use cerebro_client::terminal::ApiLoginArgs;

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
    /// SeaweedFS master node address
    #[clap(
        long, 
        short = 's', 
        default_value = "http://fs.cerebro.localhost", 
        env = "CEREBRO_FS_URL"
    )]
    pub fs_url: String,
    /// SeaweedFS master node port
    #[clap(
        long, 
        short = 'p',
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
    Ping(ApiStatusArgs),
    // Login to Cerebro API
    Login(ApiLoginArgs),
    /// Upload to SeaweedFS and registration of files with Cerebro
    Upload(UploadFileArgs),
    /// Download of files from SeaweedFS
    Download(DownloadFileArgs),
    /// Delete a file from SeaweedFS
    Delete(DeleteFileArgs),
    /// List accessible files from SeaweedFS registered with Cerebro
    List(ListFileArgs),
    /// Stage files fromn SeaweedFS periodically from stage databases registered with Cerebro
    Stage(StageFileArgs),
}

#[derive(Debug, Args)]
pub struct StageFileArgs {
    /// Stage file directory
    #[clap(long, short = 'o')]
    pub outdir: PathBuf,
}


#[derive(Debug, Args)]
pub struct DeleteFileArgs {
    /// Files identifiers to delete (CerebroFS)
    #[clap(long, short = 'f', num_args(0..))]
    pub file_ids: Vec<String>,
    /// Team name for file query (CerebroAPI)
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Database name for file query (CerebroAPI)
    #[clap(long, short = 'd')]
    pub db_name: String,
    /// Sequence run identifier
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
}


#[derive(Debug, Args)]
pub struct UploadFileArgs {
    /// Files to register
    #[clap(long, short = 'f', num_args(0..))]
    pub files: Vec<PathBuf>,
    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: String,
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
    // Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: String,
    /// Sequence run identifier
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
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
