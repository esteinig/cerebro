use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

use crate::stack::deploy::StackConfigTemplate;

/// Cerebro: metagenomic diagnostic for clinical production
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro", version)]
pub struct App {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[clap(subcommand)]
    /// Stack configuration and deployment
    Stack(StackCommands),
}


/*
==================================
STACK CONFIGURATION AND DEPLOYMENT
==================================
*/

#[derive(Debug, Subcommand)]
pub enum StackCommands {
    /// Deploy a stack configuration from file
    Deploy(StackDeployArgs),
    /// Hash a password using the stack hash function
    Hash(StackHashPasswordArgs)
}   

#[derive(Debug, Args)]
pub struct StackDeployArgs {

    /// Stack configuration name
    #[clap(long, short = 'n', env = "CEREBRO_STACK_CONFIG_NAME")]
    pub name: String,
    /// Configured stack output directory
    #[clap(long, short = 'o', env = "CEREBRO_STACK_CONFIG_DIR")]
    pub outdir: PathBuf,
    /// Stack configuration template
    #[clap(long, short = 'c', env = "CEREBRO_STACK_CONFIG_TEMPLATE")]
    pub config: Option<StackConfigTemplate>,
    /// Stack configuration file (.toml)
    #[clap(long, short = 'f', env = "CEREBRO_STACK_CONFIG_FILE")]
    pub config_file: Option<PathBuf>,
    /// Deploy a configuration from template with interactive 
    /// prompts to fill in required arguments like admin username,
    /// email and passwords if not provided on the command-line
    #[clap(long, short = 'i')]
    pub interactive: bool,
    /// Deploy for local development with hot-reloads (unsafe in production)
    #[clap(long, short = 'd')]
    pub dev: bool,
    /// Create a `.trigger` file in the deployed repository for manual rebuilds during local development
    /// 
    /// Can be used to modify the trigger file (in the base repository) to manually
    /// launch crate compilation rather than watching for any change in `src` or 
    /// `Cargo.toml` for hot-rebuilds during development
    #[clap(long, short = 't')]
    pub trigger: bool,
    /// Clone the specific branch for this repository
    #[clap(long, short = 'b')]
    pub branch: Option<String>,
    /// Checkout repository at the commit or tag provided, overwrites branch if provided
    #[clap(long, short = 'r')]
    pub revision: Option<String>,
    /// Change the subdomain for multiple concurrent stack deployments 
    /// 
    /// This will configure Traefik to deploy to {app,api}.{subdomain}.{domain}
    /// so that multiple subdomains can be configured on the fly during deployment
    /// for example to `app.dev.cerebro.localhost` or `app.demo.cerebro.localhost`
    #[clap(long, short = 's')]
    pub subdomain: Option<String>,
    /// Public or SSH-like repository URL for cloning into deployment
    #[clap(long, short = 'u', env = "CEREBRO_STACK_GIT_REPO_URL", default_value="git@github.com:esteinig/cerebro.git")]
    pub git_url: String,
    /// Primary file system path if using Cerebro FS
    #[clap(long)]
    pub fs_primary: Option<PathBuf>,
    /// Secondary file system path if using Cerebro FS
    #[clap(long)]
    pub fs_secondary: Option<PathBuf>,
    /// Database root username
    #[clap(long)]
    pub db_root_username: Option<String>,
    /// Database root password
    #[clap(long)]
    pub db_root_password: Option<String>,
    /// Database admin username
    #[clap(long)]
    pub db_admin_username: Option<String>,
    /// Database admin password
    #[clap(long)]
    pub db_admin_password: Option<String>,
    /// Cerebro admin email for admin login
    #[clap(long)]
    pub cerebro_admin_email: Option<String>,
    /// Cerebro admin password for admin login
    #[clap(long)]
    pub cerebro_admin_password: Option<String>,
    /// Cerebro admin full name for admin profile
    #[clap(long)]
    pub cerebro_admin_name: Option<String>,
    /// Email name used in domain registration used to configure Traefik
    #[clap(long)]
    pub traefik_domain: Option<String>,
    /// User name for Traefik web dashboard (BasicAuth)
    #[clap(long)]
    pub traefik_username: Option<String>,
    /// User password for Traefik web dashboard (BasicAuth)
    #[clap(long)]
    pub traefik_password: Option<String>,
}

#[derive(Debug, Args)]
pub struct StackHashPasswordArgs {
    /// Password to hash
    #[clap(long, short = 'p', env = "CEREBRO_STACK_HASH_PASSWORD")]
    pub password: String
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
