use std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

/// Cerebro: metagenomic diagnostic for clinical production
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro", version)]
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
    #[clap(subcommand)]
    /// Pipeline processing
    Pipeline(PipelineCommands),
    #[clap(subcommand)]
    /// Analytical routines
    Analysis(AnalysisCommands),
    #[clap(subcommand)]
    /// Report configurations and compiler
    Report(ReportCommands),
    #[clap(subcommand)]
    /// Utility tools for pipeline and stack
    Tools(ToolCommands),
    #[clap(subcommand)]
    /// Stack configuration and deployment
    Stack(StackCommands),
    #[clap(subcommand)]
    /// Application programming interface client
    Api(ApiCommands),
}



#[derive(Debug, Subcommand)]
pub enum ToolCommands {
    #[clap(subcommand)]
    /// Utility tools
    Utils(UtilCommands),
    #[clap(subcommand)]
    /// UMI processing tools
    Umi(UmiCommands),
    #[clap(subcommand)]
    /// Virus processing tools
    Virus(VirusCommands),

}   

#[derive(Debug, Subcommand)]
pub enum UtilCommands {
    /// Subset a database sequence file
    Subset(UtilsSubsetArgs),
    /// Anonymize read headers
    Anonymize(UtilsAnonymizeArgs),
    /// Split a FASTA file into sequence FASTA
    Split(UtilsSplitArgs),
}

#[derive(Debug, Args)]
pub struct UtilsSubsetArgs {
    /// Input database file for subset (.fasta)
    #[clap(long, short = 'f')]
    pub fasta: PathBuf,
    /// Input Mash output file 
    #[clap(long, short = 'm')]
    pub mash: PathBuf,
    /// Minimum hash identity
    #[clap(long, short = 'i')]
    pub min_identity: f64,
    /// Minimum shared hashes
    #[clap(long, short = 's')]
    pub min_shared_hashes: u32,
    /// Group subset sequences by field index in seq id
    /// 
    /// Pick the best performing from preselection after
    /// grouping. For example, sequences can be grouped
    /// by taxonomic identifier using the field index
    /// (0-based) after splitting sequence identifier with
    /// the group separator. In a subset FASTA / Mash screen
    /// result file where sequences have identifiers such as
    /// `kraken:taxid|561307|NZ_CP076105.1`
    ///  a grouping using the tax identifier can be done:
    /// --group-sep "|" and --group-index 1
    #[clap(long, short = 'g')]
    pub group_index: Option<usize>,
    /// Group subset sequences by field specifier in description
    ///
    /// Pick the best performing from preselection after
    /// grouping. For example, sequences can be grouped
    /// by taxonomic identifier using the field specifier
    /// (starts with) after splitting sequence identifier with
    /// the group separator. In a subset FASTA / Mash screen
    /// result file where sequences have desriptions such as
    /// `>seq_id taxid=561307 assembly=NZ_CP076105.1`
    /// a grouping using the tax identifier can be done:
    /// --group-sep " " and --group-by "taxid="
    #[clap(long, short = 'b')]
    pub group_by: Option<String>,
    /// Group seperator for subset on gorup index
    #[clap(long, short = 'd', default_value=" ")]
    pub group_sep: String,
    /// Output deduplicated read files (length and order must match input files) (.fastq)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
}


#[derive(Debug, Args)]
pub struct UtilsAnonymizeArgs {
    /// Input filepath(s) (fa, fq, gz, bz).
    ///
    /// For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two
    /// files consecutively `-i r1.fq r2.fq`. NOTE: Read identifiers for paired-end Illumina reads
    /// are assumed to be the same in forward and reverse read files (modern format) without trailing
    /// read orientations e.g. `/1` and `/2`. If you are using legacy identifiers, reads in the anonymized
    /// output may be unpaired.
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output filepath(s) with anonymized read identifiers.
    ///
    /// For paired Illumina you may either pass this flag twice `-o r1.fq -o r2.fq` or give two
    /// files consecutively `-o r1.fq r2.fq`. NOTE: The order of the pairs is assumed to be the
    /// same as that given for --input.
    #[clap(long, short = 'o', num_args(0..))]
    pub output: Vec<PathBuf>,
    /// Illumina fake header
    #[clap(long, short = 'f')]
    pub illumina: bool,
}

#[derive(Debug, Args)]
pub struct UtilsSplitArgs {
    /// Input filepath(s) (fa, fq, gz, bz).
    ///
    /// For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two
    /// files consecutively `-i r1.fq r2.fq`. NOTE: Read identifiers for paired-end Illumina reads
    /// are assumed to be the same in forward and reverse read files (modern format) without trailing
    /// read orientations e.g. `/1` and `/2`. If you are using legacy identifiers, reads in the anonymized
    /// output may be unpaired.
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output filepath(s) with anonymized read identifiers.
    ///
    /// For paired Illumina you may either pass this flag twice `-o r1.fq -o r2.fq` or give two
    /// files consecutively `-o r1.fq r2.fq`. NOTE: The order of the pairs is assumed to be the
    /// same as that given for --input.
    #[clap(long, short = 'o')]
    pub outdir: PathBuf,
}

#[derive(Debug, Subcommand)]
pub enum UmiCommands {
    /// Trim sample index from start of the molecular identifier sequence
    Trim(UmiTrimArgs),
    /// Check if reads contain the molecular identifier sequence indicated in header
    Check(UmiCheckArgs),
    /// Prepend molecular identifier to matching sequence read
    PrepareCalib(UmiPrepareCalibArgs),
    /// Check the number of naive clusters of molecular identifiers
    DedupNaive(UmiDedupNaiveArgs),
    /// Check and summarize a Calib cluster report 
    DedupCalib(UmiDedupCalibArgs),
}


#[derive(Debug, Args)]
pub struct UmiDedupCalibArgs {
    /// Input read files for deduplication (.fastq)
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output deduplicated read files  (length and order must match input files) (.fastq)
    #[clap(long, short = 'o', num_args(0..))]
    pub output: Vec<PathBuf>,
    /// Input cluster report from Calib output 
    #[clap(long, short = 'c')]
    pub clusters: PathBuf,
    /// Output cluster summary file (.tsv)
    #[clap(long, short = 's')]
    pub summary: Option<PathBuf>,
    /// Output deduplicated read identifier file (.tsv)
    #[clap(long, short = 'r')]
    pub identifiers: Option<PathBuf>,
}


#[derive(Debug, Args)]
pub struct UmiCheckArgs {
    /// Input read files (.fastq)
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output summary file (.tsv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
}

#[derive(Debug, Args)]
pub struct UmiDedupNaiveArgs {
    /// Input read files for deduplication (.fastq)
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Output deduplicated read files (length and order must match input files) (.fastq)
    #[clap(long, short = 'o', num_args(0..))]
    pub output: Vec<PathBuf>,
    /// Input reads with UMI in last field of read identifier (:) (.fastq)
    /// 
    /// Deduplicates reads with UMI prepended: clusters identical reads, and selects
    /// the best average read quality read (identfier) from each cluster to be subset 
    /// from the input files. Usually the forward reads of input files should be used.
    #[clap(long, short = 'f')]
    pub reads: PathBuf,
    /// Output cluster summary (.tsv)
    #[clap(long, short = 's')]
    pub summary: Option<PathBuf>,
    /// Do not extract UMI sequence from read identifier and prepend
    /// 
    /// If activated, only the read sequence is used for deduplication.
    #[clap(long, short = 'n')]
    pub no_umi: bool,
}

#[derive(Debug, Args)]
pub struct UmiTrimArgs {
    /// INDEX + UMI read file input (.fastq)
    #[clap(long, short = 'u')]
    pub umi: PathBuf,
    /// UMI read file output (.fastq)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Index length to trim from sequence
    #[clap(long, short = 't', default_value="8")]
    pub trim: usize,
}

#[derive(Debug, Args)]
pub struct UmiPrepareCalibArgs {
    /// Input read file to prepend UMI from read identifier (.fastq)
    #[clap(long, short = 'i')]
    pub input: PathBuf,
    /// Output reads with prepended UMI (.fastq)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
}


#[derive(Debug, Subcommand)]
pub enum VirusCommands {
    /// Process the outputs of the scan-remap-consensus alignment module
    ScanRemap(VirusScanRemapArgs)
}


#[derive(Debug, Args)]
pub struct VirusScanRemapArgs {
    /// Sample identifier
    #[clap(long, short = 'i')]
    pub id: String,
    /// Database name
    #[clap(long, short = 'd')]
    pub db: String,
    /// Vircov result file from the scanning module
    #[clap(long, short = 's')]
    pub scan: PathBuf,
    /// Vircov result file (concatenated) from the remap module
    #[clap(long, short = 'r')]
    pub remap: PathBuf,
    /// Consensus sequences (concatenated) from the consensus module
    #[clap(long, short = 'c')]
    pub consensus: Option<PathBuf>,
    /// Samtools coverage data (concatenated) from the remapping module
    #[clap(long, short = 'c')]
    pub coverage: Option<PathBuf>,
    /// Annotation options preset
    #[clap(long, short = 'a')]
    pub annotation: Option<String>,
    /// Output summary table (.tsv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Header for output table
    #[clap(long, short = 'H')]
    pub header: bool,
}

/*
==================================
STACK CONFIGURATION AND DEPLOYMENT
==================================
*/

#[derive(Debug, Subcommand)]
pub enum StackCommands {
    /// Deploy a stack configuration
    Deploy(StackDeployArgs),
    /// Run the stack server
    RunServer(StackServerArgs)
}   

#[derive(Debug, Args)]
pub struct StackDeployArgs {
    /// Stack configuration file (.toml)
    #[clap(long, short = 'c', env = "STACK_CONFIG_FILE")]
    pub config: PathBuf,
    /// Cerebro repository directory for development stack 
    #[clap(long, short = 'd')]
    pub dev: Option<PathBuf>,
    /// Configured stack output directory
    #[clap(long, short = 'o', env = "STACK_CONFIG_DIR")]
    pub outdir: PathBuf,
}


#[derive(Debug, Args)]
pub struct StackServerArgs {

    /// Configuration file or environmental variable
    #[clap(long, short = 'c', env = "CEREBRO_CONFIG_FILE")]
    pub config: PathBuf,
    /// Host address
    #[clap(long, short = 'H', default_value="127.0.0.1", env = "CEREBRO_SERVER_HOST")]
    pub host: String,
    /// Host port
    #[clap(long, short = 'P', default_value="8080", env = "CEREBRO_SERVER_PORT")]
    pub port: u16,
    /// Threads
    #[clap(long, short = 't', default_value="4")]
    pub threads: usize,
}


/*
===============================
PIPELINE PARSERS AND PROCESSORS
===============================
*/

#[derive(Debug, Subcommand)]
pub enum PipelineCommands {
    /// Parse and process pipeline results
    Process(PipelineProcessArgs),
    /// Quality control table from processed results
    Quality(PipelineQualityArgs),
    /// Create a sample sheet from the input directory
    SampleSheet(PipelineSampleSheetArgs)
}

#[derive(Debug, Args)]
pub struct PipelineProcessArgs {
    
    /// Pipeline sample results directory for processing
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Taxonomy directory containing 'names.dmp' and 'nodes.dmp' used for classification (NCBI)
    /// 
    /// Must be present for classification processing, otherwise only
    /// the quality control module is proccessed.
    #[clap(long, short = 't')]
    pub taxonomy: Option<PathBuf>,

    /// Output file of processed sample database model (.json)
    #[clap(long, short = 'o')]
    pub output: PathBuf,

    /// Optional sample identifier to use for output instead of result directory name
    #[clap(long, short = 's')]
    pub sample_id: Option<String>,

}


#[derive(Debug, Args)]
pub struct PipelineQualityArgs {
    
    /// Processed pipeline result samples (.json)
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,

    /// Output file for quality control table (.tsv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,

    /// Include header in quality control table
    #[clap(long, short = 'H')]
    pub header: bool,
    
    /// Input mass of ERCC or EDCC for biomass calculations (in pg)
    #[clap(long, short = 'e')]
    pub ercc_mass: Option<f64>,

}


#[derive(Debug, Args)]
pub struct PipelineSampleSheetArgs {
    
    /// Processed pipeline result samples (.json)
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,

    /// Output sample sheet file (.csv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,

    /// Sample glob - pattern matching to find paired samples
    /// 
    /// The glob string to specify paired sample extensions should be in format: {forward_extension,reverse_extension}
    /// where the wildcard specifies the sample identifier, for example: *_{R1,R2}.fastq.gz will match files
    /// "Sample1_R1.fastq.gz" and "Sample1_R2.fastq.gz" to sample identifier "Sample1"
    #[clap(long, short = 'g', default_value = "*_{R1_001,R2_001}.fastq.gz")]
    pub glob: String,

    /// Run identifier - if not provided uses input directory name
    /// 
    /// If you want to fill in custom run identifiers for each sample, 
    /// you can provide an empty string ("") and edit the sample sheet 
    /// manually.
    #[clap(long, short = 'r')]
    pub run_id: Option<String>,

    /// Run date - if not provided uses current date (YYYYMMDD)
    /// 
    /// If you want to fill in custom run dates for each sample, 
    /// you can provide an empty string ("") and edit the sample sheet 
    /// manually.
    #[clap(long, short = 'd')]
    pub run_date: Option<String>,

    /// Sample group - if not provided sample group designation is an empty string
    /// 
    /// Sample groups can be specified manually for larger runs containing
    /// sampels from multiple experimental groups - these are later available
    /// in the front-end application
    #[clap(long, short = 's')]
    pub sample_group: Option<String>,

    /// Sample type - if not provided sample type designation is an empty string
    /// 
    /// Sample types can be specified manually for larger runs containing
    /// sampels from multiple biological sources- these are later available
    /// in the front-end application
    #[clap(long, short = 't')]
    pub sample_type: Option<String>,

    /// ERCC input mass in picogram - if not provided input mass is 0
    /// 
    /// In the validation experiments, we test different input masses per sample.
    /// Generally not needed and can be overwritten with options in the fixed
    /// workflow settings later. Set to 25 pg for standard ERCC.
    #[clap(long, short = 'e')]
    pub ercc_input: Option<f64>,

    /// Allow symlink target reading for glob file walking
    #[clap(long, short = 'l')]
    pub symlinks: bool

}

/*
========================================
APPLICATION PROGRAMMING INTERFACE CLIENT
========================================
*/

#[derive(Debug, Subcommand)]
pub enum ApiCommands {
    /// Login user for authentication token 
    Login(ApiLoginArgs),
    /// Ping the server as authenticated user 
    Ping(ApiPingArgs),
    /// Ping the server as unauthenticated user
    Status(ApiStatusArgs),
    /// Upload pipeline outputs to database
    Upload(ApiUploadArgs),
    /// Summary of taxa evidence for requested models
    Taxa(ApiTaxaArgs),

    /// CRUD operations for teams
    #[clap(subcommand)]
    Team(ApiTeamCommands),
    #[clap(subcommand)]
    /// CRUD operations for team databases
    Database(ApiDatabaseCommands),
    #[clap(subcommand)]
    /// CRUD operations for team database projects
    Project(ApiProjectCommands),
}

#[derive(Debug, Args)]
pub struct ApiPingArgs {
}

#[derive(Debug, Args)]
pub struct ApiStatusArgs {
}

#[derive(Debug, Args)]
pub struct ApiLoginArgs {
    /// Registered user email
    #[clap(long, short = 'e', env = "CEREBRO_USER_EMAIL")]
    pub email: String,
    /// Registered user password
    #[clap(long, short = 'p', env = "CEREBRO_USER_PASSWORD")]
    pub password: Option<String>,
}


#[derive(Debug, Args)]
pub struct ApiUploadArgs {

    /// Processed pipeline sample model (.json)
    #[clap(long, short = 'i')]
    pub input: PathBuf,

    /// Pipeline sample sheet (.csv)
    #[clap(long, short = 's')]
    pub sample_sheet: PathBuf,
    
    /// Pipeline configuration (.json)
    #[clap(long, short = 'c')]
    pub pipeline_config: PathBuf,

    /// Team name for model upload
    #[clap(long, short = 't')]
    pub team_name: String,

    /// Project name for model upload
    #[clap(long, short = 'p')]
    pub project_name: String,

    /// Database name for model upload
    #[clap(long, short = 'd')]
    pub db_name: Option<String>,

    /// Output model as file (.json)
    #[clap(long, short = 'o')]
    pub output: Option<PathBuf>,
}


#[derive(Debug, Args)]
pub struct ApiTaxaArgs {

    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,

    /// Project name for model query
    #[clap(long, short = 'p')]
    pub project_name: String,

    /// Output summary (.csv)
    #[clap(long, short = 'o')]
    pub output: PathBuf,

    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: Option<String>,

    /// Filter JSON that satisfies `CerebroFilterConfig` schema
    #[clap(long, short = 'f')]
    pub filter_config: Option<PathBuf>,

    /// Run identifiers to filter results
    #[clap(long, short = 'r', num_args(0..))]
    pub run_ids: Option<Vec<String>>,

    /// Sample identifiers to filter results
    #[clap(long, short = 's', num_args(0..))]
    pub sample_ids: Option<Vec<String>>,

    /// Workflow identifiers to filter results
    #[clap(long, short = 'w', num_args(0..))]
    pub workflow_ids: Option<Vec<String>>,

    /// Workflow mnemnonic names to filter results
    #[clap(long, short = 'n', num_args(0..))]
    pub workflow_names: Option<Vec<String>>,
}

#[derive(Debug, Args)]
pub struct ApiTeamArgs {

    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,

    /// Project name for model query
    #[clap(long, short = 'p')]
    pub team_descriptions: String,

}

#[derive(Debug, Subcommand)]
pub enum ApiTeamCommands {
    // Create a new team
    Create(ApiTeamCreateArgs)
}


#[derive(Debug, Subcommand)]
pub enum ApiDatabaseCommands {
    // Create a new team
    Create(ApiDatabaseCreateArgs)
}


#[derive(Debug, Subcommand)]
pub enum ApiProjectCommands {
    // Create a new team
    Create(ApiProjectCreateArgs)
}


#[derive(Debug, Args)]
pub struct ApiTeamCreateArgs {

    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,

    /// Project name for model query
    #[clap(long, short = 'p')]
    pub team_descriptions: String,

}

#[derive(Debug, Args)]
pub struct ApiDatabaseCreateArgs {

    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,

    /// Project name for model query
    #[clap(long, short = 'p')]
    pub project_name: String,

    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: Option<String>,

}

#[derive(Debug, Args)]
pub struct ApiProjectCreateArgs {

    /// Team name for model query
    #[clap(long, short = 't')]
    pub team_name: String,

    /// Database name for model query
    #[clap(long, short = 'd')]
    pub db_name: String,

    /// Project name for model query
    #[clap(long, short = 'n')]
    pub project_name: String,

    /// Project name for model query
    #[clap(long, short = 'i')]
    pub project_description: String,

}



/*
=======================
REPORT TEMPLATE ENGINE
======================
*/

#[derive(Debug, Subcommand)]
pub enum ReportCommands {
    /// Generate clinical reports from data and templates
    Compile(ReportCompileArgs),
}

#[derive(Debug, Args)]
pub struct ReportCompileArgs {
    /// Output report LaTeX template
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Compile LaTeX template into output report
    #[cfg(feature = "pdf")]
    #[clap(long, short = 'p')]
    pub pdf: bool,
    /// Base template configuration file (.toml)
    #[clap(long, short = 'c')]
    pub base_config: PathBuf,
    /// Sample template configuration file (.toml)
    #[clap(long, short = 's')]
    pub sample_config: Option<Vec<PathBuf>>,
    /// Complete: patient template configuration file (.toml)
    #[clap(long, short = 'P')]
    pub patient_config: Option<PathBuf>,
    /// Partial: patient header template configuration file (.toml)
    #[clap(long, short = 'H')]
    pub patient_header_config: Option<PathBuf>,
    /// Partial: patient result template configuration file (.toml)
    #[clap(long, short = 'R')]
    pub patient_result_config: Option<PathBuf>,
}



#[derive(Debug, Subcommand)]
pub enum AnalysisCommands {

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
