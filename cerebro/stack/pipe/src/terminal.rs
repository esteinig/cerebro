use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

/// Cerebro: production stack server
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-pipeline", version)]
pub struct App {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
      

    #[clap(subcommand)]
    /// Parse and process pipeline results
    Process(ProcessCommands),

    #[clap(subcommand)]
    /// Pipeline output summary tables
    Tables(TablesCommands),

    #[clap(subcommand)]
    /// Internal pipeline utilities and sample sheets
    Tools(ToolsCommands),
}




#[derive(Debug, Subcommand)]
pub enum SheetCommands {
}


#[derive(Debug, Subcommand)]
pub enum TablesCommands {
    /// Quality control tables
    QualityControl(QualityControlTableArgs),
}

#[derive(Debug, Args)]
pub struct QualityControlTableArgs {
    /// Quality control summaries (.json)
    #[clap(long, short = 'j', num_args(0..))]
    pub json: Vec<PathBuf>,
    /// Summary read quality and alignment table
    #[clap(long, short = 'r')]
    pub reads: PathBuf,
    /// Controls alignment table
    #[clap(long, short = 'c')]
    pub controls: PathBuf,
    /// Background alignment table
    #[clap(long, short = 'b')]
    pub background: PathBuf,
}

#[derive(Debug, Subcommand)]
pub enum ToolsCommands {
    /// Scan read files and compute summary metrics
    ScanReads(ScanReadArgs),
    /// Create a sample-sheet for read input to pipelines
    SampleSheet(SampleSheetArgs),
    /// Check if a file contains UMIs in read identifier
    UmiCheck(UmiCheckArgs),
    /// Trim the I2 files from UMI demultiplexed Illumina data
    UmiTrim(UmiTrimArgs),
    /// Prepare Calib UMI deduplication input files
    UmiPrepCalib(UmiPrepareCalibArgs),
    /// Deduplicate UMI barcodes with Calib
    UmiDedupCalib(UmiDedupCalibArgs),
    /// Deduplicate UMI barcodes using a naive clustering algorithm
    UmiDedupNaive(UmiDedupNaiveArgs),
}

#[derive(Debug, Args)]
pub struct ScanReadArgs {
    /// Input directory containing the output files for the panviral pipeline
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Sample identifier to parse, directory basename by default
    #[clap(long, short = 'j')]
    pub json: Option<PathBuf>,
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
    /// Use only the first x base pairs of a read for deduplication or entire read if shorter
    #[clap(long, short = 'x')]
    pub head: Option<usize>,
    /// Output deduplication report (.json)
    #[clap(long, short = 'j')]
    pub json: Option<PathBuf>,
    /// Output cluster summary (.tsv)
    #[clap(long, short = 'c')]
    pub clusters: Option<PathBuf>,
    /// Delimiter in tail of read ids to extract UMI sequence
    #[clap(long, short = 'd', default_value=":")]
    pub delim: String,
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


#[derive(Debug, Args)]
pub struct SampleSheetArgs {
    /// Input directory to collect sample files from
    #[clap(long, short = 'd', num_args(0..))]
    pub directory: Vec<PathBuf>,
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
    /// Sample prefix file - only collect samples starting with prefix
    /// 
    /// One sample prefix per line to collect using the directory and glob search.
    #[clap(long, short = 'p')]
    pub prefix: Option<PathBuf>,
    /// Run identifier - if not provided uses input directory name
    /// 
    /// If you want to fill in custom run identifiers for each sample, 
    /// you can provide an empty string ("") and edit the sample sheet 
    /// manually.
    #[clap(long, short = 'i')]
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
    /// workflow settings later. Set to 25 pg for standard ERCC inputs.
    #[clap(long, short = 'e')]
    pub ercc_input: Option<f64>,
    /// Allow symlink target reading for glob file walking
    #[clap(long, short = 'l')]
    pub symlinks: bool,
    /// Find files recursively in input path
    #[clap(long, short = 'r')]
    pub recursive: bool,
    /// Allow non-unique sample identifiers in sample sheet (caution)
    #[clap(long, short = 'u')]
    pub not_unique: bool,
    /// Single reads for Oxford Nanopore Technology 
    #[clap(long)]
    pub ont: bool
}


#[derive(Debug, Subcommand)]
pub enum ProcessCommands {
    /// Process panviral enrichment outputs
    Panviral(ProcessArgs),
    /// Process pathogen detection outputs
    Pathogen(ProcessArgs),
    /// Process quality control outputs
    Quality(ProcessArgs),
}


#[derive(Debug, Args)]
pub struct ProcessArgs {
    /// Input directory containing the output files for the pathogen pipeline
    #[clap(long, short = 'i', default_value=".")]
    pub input: PathBuf,
    /// Sample identifier to parse, directory basename by default
    #[clap(long, short = 's')]
    pub id: Option<String>,
    /// Output file of processed quality control data
    #[clap(long, short = 'q')]
    pub qc: PathBuf,
    /// Parse the background alignment from the quality control module variant
    /// 
    /// Combined reference for detecting organism, synthetic, internal controls 
    /// and other relevant background and deplete in a single alignment before
    /// other quality control steps
    #[clap(long, short = 'a')]
    pub background: bool,

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