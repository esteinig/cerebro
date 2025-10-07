use std::path::PathBuf;
use cerebro_model::api::{watchers::model::WatcherFormat, towers::model::Pipeline};
use clap::{ArgGroup, Args, Parser, Subcommand};

/// Cerebro: production stack server
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-client", version)]
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
    /// User team name or identifier for requests that require team data acess 
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
    /// API token file - can be set from environment variable
    #[clap(
        long, 
        short = 'f', 
        env = "CEREBRO_API_TOKEN_FILE"
    )]
    pub token_file: Option<PathBuf>,
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
    
    /// Login user for authentication token 
    Login(LoginArgs),
    /// Ping the server as authenticated user 
    PingServer(PingArgs),
    /// Ping the server as unauthenticated user
    PingStatus(StatusArgs),

    /// Process pathogen-detection pipeline outputs for upload to database
    CreatePathogen(CreatePathogenArgs),
    /// Process panviral-enrichment pipeline outputs for upload to database
    CreatePanviral(CreatePanviralArgs),

    /// Upload pipeline outputs to database (production from staged files)
    UploadPathogen(UploadPathogenArgs),
    /// Upload pipeline outputs to database (production from staged files)
    UploadPanviral(UploadPanviralArgs),
    
    /// Upload processed model to database
    UploadModels(UploadModelsArgs),
    /// Upload processed model to database
    DownloadModels(DownloadModelsArgs),

    /// Summary of quality control for requested models or samples
    GetQualityControlTable(QualityControlTableArgs),
    /// Summary of pathogen detection for requested models or samples
    GetPathogenDetectionTable(PathogenDetectionTableArgs),
    
    /// Taxon history query for data and host regression
    GetTaxonHistory(TaxonHistoryArgs),


    /// CRUD operations for training collections
    #[clap(subcommand)]
    Training(TrainingCommands),

    /// CRUD operations for production watchers
    #[clap(subcommand)]
    Watcher(WatcherCommands),
    /// CRUD operations for production pipelines
    #[clap(subcommand)]
    Tower(TowerCommands),
    /// CRUD operations for staging production files
    #[clap(subcommand)]
    Stage(StageCommands),

    /// CRUD operations for teams
    #[clap(subcommand)]
    Team(TeamCommands),
    #[clap(subcommand)]
    /// CRUD operations for team databases
    Database(DatabaseCommands),
    #[clap(subcommand)]
    /// CRUD operations for team database projects
    Project(ProjectCommands),

    /// Faktory job scheduler
    #[clap(subcommand)]
    Jobs(JobsCommands),
}

#[derive(Debug, Args)]
pub struct PingArgs {
}

#[derive(Debug, Args)]
pub struct StatusArgs {
}

#[derive(Debug, Args)]
pub struct LoginArgs {
    /// Registered user email
    #[clap(long, short = 'e', env = "CEREBRO_USER_EMAIL")]
    pub email: String,
    /// Registered user password
    #[clap(long, short = 'p', env = "CEREBRO_USER_PASSWORD")]
    pub password: Option<String>,
}


#[derive(Debug, Args)]
pub struct UploadPathogenArgs {
    /// Processed pipeline quality control module (.json)
    #[clap(long)]
    pub quality: PathBuf,
    /// Processed pipeline pathogen detection module (.json)
    #[clap(long)]
    pub pathogen: PathBuf,
    /// Taxonomy directory containing 'nodes.dmp' and 'names.dmp'
    #[clap(long)]
    pub taxonomy: PathBuf,
    /// Raise error if taxid was not found in taxonomy
    #[clap(long)]
    pub strict: bool,
    /// GTDB - collapse paraphyletic species into the same species - may reassign taxonomic identifiers 
    #[clap(long)]
    pub gtdb_break_monophyly: bool,
    /// Run identifier if sample sheet is not provided
    #[clap(long)]
    pub run_id: Option<String>,
    /// Run date (YYYYMMDD) if sample sheet is not provided
    #[clap(long)]
    pub run_date: Option<String>,
    /// Pipeline sample sheet (.csv)
    #[clap(long)]
    pub sample_sheet: Option<PathBuf>,
    /// Pipeline configuration (.json)
    #[clap(long)]
    pub pipeline_config: Option<PathBuf>,
    /// Staged file metadata from production pipeline (.json)
    #[clap(long, short = 'j')]
    pub stage_json: Option<PathBuf>,
    /// Output database model as file (.json)
    #[clap(long, short = 'o')]
    pub model: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct CreatePathogenArgs {
    /// Output database model as file (.json)
    #[clap(long, short = 'o')]
    pub model: PathBuf,
    /// Processed pipeline quality control module (.json)
    #[clap(long)]
    pub quality: PathBuf,
    /// Processed pipeline pathogen detection module (.json)
    #[clap(long)]
    pub pathogen: PathBuf,
    /// Taxonomy directory containing 'nodes.dmp' and 'names.dmp'
    #[clap(long)]
    pub taxonomy: PathBuf,
    /// Raise error if taxid was not found in taxonomy
    #[clap(long)]
    pub strict: bool,
    /// GTDB - collapse paraphyletic species into the same species - may reassign taxonomic identifiers 
    #[clap(long)]
    pub gtdb_break_monophyly: bool,
    /// Run identifier if sample sheet is not provided
    /// otherwise defaults to placeholder for now
    #[clap(long)]
    pub run_id: Option<String>,
    /// Run date (YYYYMMDD) if sample sheet is not provided
    /// otherwise defaults to current date
    #[clap(long)]
    pub run_date: Option<String>,
    /// Pipeline sample sheet (.csv)
    #[clap(long)]
    pub sample_sheet: Option<PathBuf>,
    /// Pipeline configuration (.json)
    #[clap(long)]
    pub pipeline_config: Option<PathBuf>,
}



#[derive(Debug, Args)]
pub struct UploadPanviralArgs {
    /// Processed pipeline quality control module (.json)
    #[clap(long)]
    pub quality: PathBuf,
    /// Processed pipeline panviral detection module (.json)
    #[clap(long)]
    pub panviral: PathBuf,
    /// Taxonomy directory containing 'nodes.dmp' and 'names.dmp'
    #[clap(long)]
    pub taxonomy: PathBuf,
    /// Raise error if taxid was not found in taxonomy
    #[clap(long)]
    pub strict: bool,
    /// GTDB - collapse paraphyletic species into the same species - may reassign taxonomic identifiers 
    #[clap(long)]
    pub gtdb_break_monophyly: bool,
    /// Run identifier if sample sheet is not provided
    /// otherwise defaults to placeholder for now
    #[clap(long)]
    pub run_id: Option<String>,
    /// Run date (YYYYMMDD) if sample sheet is not provided
    /// otherwise defaults to current date
    #[clap(long)]
    pub run_date: Option<String>,
    /// Pipeline sample sheet (.csv)
    #[clap(long)]
    pub sample_sheet: Option<PathBuf>,
    /// Pipeline configuration (.json)
    #[clap(long)]
    pub pipeline_config: Option<PathBuf>,
    /// Staged file metadata from production pipeline (.json)
    #[clap(long, short = 'j')]
    pub stage_json: Option<PathBuf>,
    /// Output database model as file (.json)
    #[clap(long, short = 'o')]
    pub model: Option<PathBuf>,
}



#[derive(Debug, Args)]
pub struct CreatePanviralArgs {
    /// Output database model as file (.json)
    #[clap(long, short = 'o')]
    pub model: PathBuf,
    /// Processed pipeline quality control module (.json)
    #[clap(long)]
    pub quality: PathBuf,
    /// Processed pipeline panviral detection module (.json)
    #[clap(long)]
    pub panviral: PathBuf,
    /// Taxonomy directory containing 'nodes.dmp' and 'names.dmp'
    #[clap(long)]
    pub taxonomy: PathBuf,
    /// Raise error if taxid was not found in taxonomy
    #[clap(long)]
    pub strict: bool,
    /// GTDB - collapse paraphyletic species into the same species - may reassign taxonomic identifiers 
    #[clap(long)]
    pub gtdb_break_monophyly: bool,
    /// Run identifier if sample sheet is not provided
    /// otherwise defaults to placeholder for now
    #[clap(long)]
    pub run_id: Option<String>,
    /// Run date (YYYYMMDD) if sample sheet is not provided
    /// otherwise defaults to current date
    #[clap(long)]
    pub run_date: Option<String>,
    /// Pipeline sample sheet (.csv)
    #[clap(long)]
    pub sample_sheet: Option<PathBuf>,
    /// Pipeline configuration (.json)
    #[clap(long)]
    pub pipeline_config: Option<PathBuf>,
}


#[derive(Debug, Args)]
pub struct UploadModelsArgs {
    /// Processed database models(.json)
    #[clap(long, short = 'm', num_args(1..))]
    pub models: Vec<PathBuf>
}

#[derive(Debug, Args)]
pub struct DownloadModelsArgs {
    /// Output directory for models in the specified collection (.json)
    #[clap(long, short = 'o')]
    pub outdir: PathBuf
}

#[derive(Debug, Args)]
pub struct QualityControlTableArgs {
    /// Output summary table (.csv)
    #[clap(long, short = 'o')]
    pub output: Option<PathBuf>,
    /// Sample identifiers to limit query
    #[clap(long, short = 's', num_args(1..))]
    pub sample_ids: Option<Vec<String>>
}


#[derive(Debug, Args)]
pub struct PathogenDetectionTableArgs {
    /// Output summary table (.csv)
    #[clap(long, short = 'o')]
    pub output: Option<PathBuf>,
    /// Sample identifiers to limit query
    #[clap(long, short = 's', num_args(1..))]
    pub sample_ids: Option<Vec<String>>,
    /// Collapse taxon variants into single species (WARNING: assigns new taxid to the collapsed taxon)
    #[clap(long, short = 'c')]
    pub collapse_variants: bool
}


#[derive(Debug, Args)]
pub struct TaxonHistoryArgs {
    /// GTDB-like taxon label to query (example: "s__Staphylococcus aureus")
    #[clap(long, short = 't')]
    pub taxon_label: String,
    /// GTDB-like host label to query (example: "s__Homo sapiens")
    #[clap(long, short = 'l', default_value="s__Homo sapiens")]
    pub host_label: String,
    /// Output regression data instead of history data
    #[clap(long, short = 'r')]
    pub regression: bool,
    /// Output regression plot to this file
    #[clap(long, short = 'p')]
    pub plot: Option<PathBuf>,
}

#[derive(clap::Subcommand, Debug, Clone)]
pub enum JobsCommands {
    /// Enqueue a job to run ASAP
    LaunchJob {
        /// Worker type, e.g. "gridfs_process"
        #[arg(long)]
        kind: String,
        /// JSON string for object-form args (e.g. '{"db":"...","project":"..."}')
        #[arg(long)]
        args: String,
        /// Faktory queue (default: "default")
        #[arg(long)]
        queue: Option<String>,
        /// Retry count
        #[arg(long)]
        retry: Option<i32>,
        /// Reserve-for seconds
        #[arg(long, alias="ttr")]
        reserve_for_seconds: Option<u64>,
    },
    /// Schedule a one-shot or recurring job
    ScheduleJob {
        #[arg(long)]
        kind: String,
        /// JSON string for object-form args
        #[arg(long)]
        args: String,
        /// First run time in RFC3339 (UTC), e.g. "2025-08-25T09:00:00Z"
        #[arg(long)]
        run_at: String,
        /// Optional: recurrence period in seconds
        #[arg(long)]
        interval_seconds: Option<i64>,
        #[arg(long)]
        queue: Option<String>,
        #[arg(long)]
        retry: Option<i32>,
        #[arg(long, alias="ttr")]
        reserve_for_seconds: Option<u64>,
        #[arg(long)]
        enabled: Option<bool>,
    },

    /// Show the scheduler summary (next/previous 10)
    ScheduleSummary,

    /// Show a one-line status summary for a job id (UUID or Faktory JID)
    JobStatus {
        /// Job id: your UUID or Faktory JID
        #[arg(long)]
        id: String,
    },

    /// Print "yes"/"no" if the job is complete
    JobComplete {
        /// Job id: your UUID or Faktory JID
        #[arg(long)]
        id: String,
    },
}


#[derive(Debug, Subcommand)]
pub enum TrainingCommands {
    /// Upload and register prefetch data in the team database
    Upload(TrainingUploadArgs),
    /// List registered prefetch data by training collection
    List(TrainingListArgs),
    /// Delete registered prefetch data 
    Delete(TrainingDeleteArgs),
}


#[derive(Debug, Args)]
pub struct TrainingUploadArgs {
    /// Prefetch data to upload
    #[clap(short, long, short = 'p', num_args(1..))]
    pub prefetch: Vec<PathBuf>,
    /// Training collection
    #[clap(short, long)]
    pub collection: String,
    /// Training collection description
    #[clap(short, long)]
    pub description: String,
}

#[derive(Debug, Args)]
pub struct TrainingListArgs {
    /// Training collection
    #[clap(short, long)]
    pub collection: Option<String>,
}


#[derive(Debug, Args)]
pub struct TrainingDeleteArgs {
    /// Training collection
    #[clap(short, long)]
    pub collection: String,
    /// Training name
    #[clap(short, long)]
    pub name: Option<String>,
}

#[derive(Debug, Args)]
pub struct TeamArgs {
    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,
    /// Project name for model query
    #[clap(long, short = 'p')]
    pub team_descriptions: String,
}

#[derive(Debug, Subcommand)]
pub enum StageCommands {
    /// Register samples for processing with a production pipeline
    Register(StageRegisterArgs),
    /// List registered samples in a pipeline staging area
    List(StageListArgs),
    /// Delete registered samples from a pipeline staging area
    Delete(StageDeleteArgs),
    /// Pull staged samples from a pipeline staging area
    Pull(StagePullArgs),
}


#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
#[clap(group = ArgGroup::new("files").required(true).args(&["file_ids", "run_id"]))]
pub struct StageRegisterArgs {
    /// Registered tower identifier
    #[clap(long, short = 'i', group = "input")]
    pub tower_id: Option<String>,
    /// Pipeline registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
    /// Registered tower identifier
    #[clap(long, short = 'p')]
    pub pipeline: Pipeline,
    /// File identifiers for sample stage registration
    #[clap(long, short = 'f', group = "files")]
    pub file_id: Option<Vec<String>>,
    /// Run name for sample stage registration
    #[clap(long, short = 'r', group = "files")]
    pub run_id: Option<String>
}



#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
pub struct StagePullArgs {
    /// Registered tower identifier
    #[clap(long, short = 'i', group = "input")]
    pub tower_id: Option<String>,
    /// Pipeline registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
    /// Run name to list in staging area
    #[clap(long, short = 'o', default_value=".")]
    pub outdir: PathBuf,
    /// Delete staged samples
    #[clap(long, short = 'd')]
    pub delete: bool,
    /// Run name to list in staging area
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Sample name to list in staging area
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
pub struct StageListArgs {
    /// Registered tower identifier
    #[clap(long, short = 'i', group = "input")]
    pub tower_id: Option<String>,
    /// Pipeline registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
    /// Run name to list in staging area
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,
    /// Sample name to list in staging area
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
#[clap(group = ArgGroup::new("delete").required(true).args(&["staged_id", "run_id", "sample_id", "all"]))]
pub struct StageDeleteArgs {
    /// Registered tower identifier
    #[clap(long, short = 'i', group = "input")]
    pub tower_id: Option<String>,
    /// Pipeline registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
    /// Staged sample unique identifier to delete from staging area
    #[clap(long, short = 'u', group = "delete")]
    pub staged_id: Option<String>,
    /// Run name to delete from staging area
    #[clap(long, short = 'r', group = "delete")]
    pub run_id: Option<String>,
    /// Sample name to delete from staging area
    #[clap(long, short = 's', group = "delete")]
    pub sample_id: Option<String>,
    /// Delete all staged samples in the staging area
    #[clap(long, short = 'a', group = "delete")]
    pub all: bool,
}


#[derive(Debug, Subcommand)]
pub enum TowerCommands {
    /// Register a new production tower
    Register(TowerRegisterArgs),
    /// List registered production towers
    List(TowerListArgs),
    /// Delete a registered production tower
    Delete(TowerDeleteArgs),
    /// Ping a registered tower to update active status
    Ping(TowerPingArgs),
}

#[derive(Debug, Args)]
pub struct TowerRegisterArgs {
    /// Tower name for registration
    #[clap(long, short = 'n')]
    pub name: String,
    /// Tower location for registration
    #[clap(long, short = 'l')]
    pub location: String,
    /// Output registration for future reference (.json)
    #[clap(long, short = 'j')]
    pub json: Option<PathBuf>,
    /// Cerebro pipelines for registration, otherwise are all available
    #[clap(long, short = 'p', num_args(0..))]
    pub pipelines: Option<Vec<Pipeline>>,
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
pub struct TowerPingArgs {
    /// Tower identifier generated during registration
    #[clap(long, short = 'i', group = "input")]
    pub id: Option<String>,
    /// Tower registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>
}

#[derive(Debug, Args)]
pub struct TowerListArgs {
    /// Tower identifier generated during registration for single record listing
    #[clap(long, short = 'i')]
    pub id: Option<String>,
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json", "name", "location", "all"]))]
pub struct TowerDeleteArgs {
    /// Tower identifier generated during registration
    #[clap(long, short = 'i', group = "input")]
    pub id: Option<String>,
    /// Tower registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
    /// Tower name
    #[clap(long, short = 'r', group = "input")]
    pub name: Option<String>,
    /// Watcher location
    #[clap(long, short = 'l', group = "input")]
    pub location: Option<String>,
    /// Delete all towers (requires confirmation)
    #[clap(long, short = 'a')]
    pub all: bool,
}



#[derive(Debug, Subcommand)]
pub enum WatcherCommands {
    /// Register a new production watcher
    Register(WatcherRegisterArgs),
    /// List registered production watchers
    List(WatcherListArgs),
    /// Delete a registered production watcher
    Delete(WatcherDeleteArgs),
    /// Ping a registered watcher to update active status
    Ping(WatcherPingArgs),
}
#[derive(Debug, Args)]
pub struct WatcherRegisterArgs {
    /// Watcher name for registration
    #[clap(long, short = 'n')]
    pub name: String,
    /// Watcher location for registration
    #[clap(long, short = 'l')]
    pub location: String,
    /// Watcher input format 
    #[clap(long, short = 'f')]
    pub format: WatcherFormat,
    /// Fastq file glob when format is 'fastq' 
    #[clap(long, short = 'g')]
    pub glob: Option<String>,
    /// Output registration for future reference (.json)
    #[clap(long, short = 'j')]
    pub json: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct WatcherListArgs {
    /// Watcher identifier generated during registration
    #[clap(long, short = 'i')]
    pub id: Option<String>,
}


#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json"]))]
pub struct WatcherPingArgs {
    /// Watcher identifier generated during registration for single record listing
    #[clap(long, short = 'i', group = "input")]
    pub id: Option<String>,
    /// Watcher registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
}

#[derive(Debug, Args)]
#[clap(group = ArgGroup::new("input").required(true).args(&["id", "json", "name", "location", "all"]))]
pub struct WatcherDeleteArgs {
    /// Watcher identifier generated during registration
    #[clap(long, short = 'i', group = "input")]
    pub id: Option<String>,
    /// Watcher registration record (.json)
    #[clap(long, short = 'j', group = "input")]
    pub json: Option<PathBuf>,
    /// Watcher name
    #[clap(long, short = 'r', group = "input")]
    pub name: Option<String>,
    /// Watcher location
    #[clap(long, short = 'l', group = "input")]
    pub location: Option<String>,
    /// Delete all watchers (requires confirmation)
    #[clap(long, short = 'a')]
    pub all: bool,
}


#[derive(Debug, Subcommand)]
pub enum TeamCommands {
    // Create a new team
    Create(TeamCreateArgs)
}


#[derive(Debug, Subcommand)]
pub enum DatabaseCommands {
    // Create a new team
    Create(DatabaseCreateArgs)
}


#[derive(Debug, Subcommand)]
pub enum ProjectCommands {
    // Create a new team
    Create(ProjectCreateArgs)
}


#[derive(Debug, Args)]
pub struct TeamCreateArgs {
    /// Team name 
    #[clap(long, short = 'n')]
    pub name: String,
    /// Team description
    #[clap(long, short = 'd')]
    pub description: String,
}

#[derive(Debug, Args)]
pub struct DatabaseCreateArgs {
    /// Database name
    #[clap(long, short = 'n')]
    pub name: String,
    /// Database description
    #[clap(long, short = 'd')]
    pub description: String,
}

#[derive(Debug, Args)]
pub struct ProjectCreateArgs {
    /// Project name
    #[clap(long, short = 'n')]
    pub name: String,
    /// Project description
    #[clap(long, short = 'd')]
    pub description: String,
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
