use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

use crate::{error::WorkflowError, modules::pathogen::PathogenDetectionRank, tools::download::{Aligner, Classifier}, tools::download::{CerebroDownloader, CerebroDownloaderBuilder, CerebroIndex}};

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
    Table(TableCommands),
    #[clap(subcommand)]
    /// Internal pipeline utilities and sample sheets
    Tools(ToolsCommands),
}


#[derive(Debug, Subcommand)]
pub enum TableCommands {
    /// Quality control tables
    QualityControl(QualityControlTableArgs),
    /// Pathogen detection (profiling) tables
    PathogenDetection(PathogenDetectionTableArgs),
    /// Filter a pathogen detection table
    FilterTable(FilterTableArgs),
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


#[derive(Debug, Args)]
pub struct PathogenDetectionTableArgs {
    /// Pathogen detection summaries (.pd.json)
    #[clap(long, short = 'j', num_args(0..))]
    pub json: Vec<PathBuf>,
    /// Reference taxonomy to extract rank and name for taxonomic identifiers
    #[clap(long, short = 'x')]
    pub taxonomy: Option<PathBuf>,
    /// Output table for aggregated taxon reads per million
    #[clap(long, short = 'o')]
    pub output: PathBuf,
}


#[derive(Debug, Args)]
pub struct FilterTableArgs {
    /// Pathogen detection summaries (.json)
    #[clap(long, short = 't', num_args(0..))]
    pub table: Vec<PathBuf>,
    /// Output table for aggregated taxon reads per million
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Provide filters as JSON
    #[clap(long, short = 'f')]
    pub filter_json: Option<PathBuf>
}

#[derive(Debug, Subcommand)]
pub enum ToolsCommands {
    /// List available indices and download files for aligners and classfiers.
    Download(DownloadArgs),
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
    /// samples from multiple experimental groups - these are later available
    /// in the front-end application and therefore should have permission from
    /// governance and ethics to store such types; please check your local 
    /// regulartory framework for inclusion of this data type.
    #[clap(long, short = 'g')]
    pub sample_group: Option<String>,
    /// Sample type - if not provided sample type designation is an empty string
    /// 
    /// Sample types can be specified manually for larger runs containing
    /// samples from multiple biological sources- these are later available
    /// in the front-end application and therefore should have permission from
    /// governance and ethics to store such types; please check your local 
    /// regulartory framework for inclusion of this data type.
    #[clap(long, short = 't')]
    pub sample_type: Option<String>,
    /// Sample date - if not provided sample date designation is an empty string
    /// 
    /// Sample dates can be specified manually for larger runs containing
    /// samples from multiple biological sources- these are later available
    /// in the front-end application and therefore should have permission from
    /// governance and ethics to store such types; please check your local 
    /// regulartory framework for inclusion of this data type.
    #[clap(long, short = 's')]
    pub sample_date: Option<String>,
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


#[derive(Args, Debug, Clone)]
pub struct DownloadArgs {
    /// List available index names and exit
    #[arg(short, long)]
    pub list: bool,
    /// Index name to download
    #[arg(short, long, num_args(0..))]
    pub index: Vec<CerebroIndex>,
    /// Output directory for index download
    #[arg(short, long, default_value=".")]
    pub outdir: PathBuf,
    /// Download index for one or more aligners 
    #[arg(short, long, num_args(0..))]
    pub aligner: Option<Vec<Aligner>>,
    /// Download index for one or more classifiers 
    #[arg(short, long, num_args(0..))]
    pub classfier: Option<Vec<Classifier>>,
    /// Download timeout in minutes - increase for large files and slow connections
    #[arg(short, long, default_value="360")]
    pub timeout: u64,
    /// Download a specific version of the database, no version downloads latest
    #[arg(short, long)]
    pub version: Option<String>,
    /// Include reference sequences (.fasta) for the selected indices
    #[arg(short, long)]
    pub reference: bool,
}
impl DownloadArgs {
    /// Validates the provided arguments and builds a `CerebroDownloader` instance.
    ///
    /// This method checks the provided arguments for consistency and constructs 
    /// a `CerebroDownloader` instance based on the validated arguments.
    ///
    /// # Returns
    ///
    /// * `Result<CerebroDownloader, WorkflowError>` - Ok with the constructed CerebroDownloader instance, otherwise an error.
    /// 
    /// # Example
    ///
    /// ```
    /// use clap::Parser;
    /// 
    /// let dl_args = Download::parse();
    /// let dl = dl_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(self) -> Result<CerebroDownloader, WorkflowError> {
        
        let downloader = CerebroDownloaderBuilder::new(
            self.outdir, self.index
        )
        .classifier(self.classfier)
        .aligner(self.aligner)
        .timeout(self.timeout)
        .reference(self.reference)
        .version(self.version)
        .build()?;

        Ok(downloader)
    }
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
    /// Input directory containing the output files for the pathogen or panviral pipeline
    #[clap(long, short = 'i', default_value=".")]
    pub input: PathBuf,
    /// Input directory containing the output files for the quality control module 
    /// 
    /// If not provided will check the '--input' directory for the quality control module outputs
    #[clap(long, short = 'q', default_value=".")]
    pub input_qc: Option<PathBuf>,
    /// Sample identifier to parse, directory basename by default
    #[clap(long, short = 's')]
    pub id: Option<String>,
    /// Output file of processed quality control data
    #[clap(long, short = 'q')]
    pub qc: Option<PathBuf>,
    /// Output file of processed pathogen detection data
    #[clap(long, short = 'p')]
    pub pathogen: Option<PathBuf>,
    /// Output file of processed panviral data
    #[clap(long, short = 'v')]
    pub panviral: Option<PathBuf>,
    /// Taxonomy for metagenome assembly contig typing (BLAST) when using LCA
    #[clap(long, short = 't')]
    pub taxonomy_directory: Option<PathBuf>,
    /// Compute LCA for taxid reassignment from BLAST hits
    #[clap(long, short = 'l')]
    pub blast_lca: bool,
    /// If read scanning outputs fail to be detected (input/output scanning) include sample with zero reads
    #[clap(long, short = 'f')]
    pub qc_fail_ok: bool,
    /// Paired-end reads (same identifier for each mate) used as input
    /// 
    /// Classifiers like Kraken2 or Metabuli output reads as unique reads
    /// identifiers that were classified, therefore the read counts are
    /// usually half of what they are expected.
    #[clap(long)]
    pub paired_end: bool,

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
