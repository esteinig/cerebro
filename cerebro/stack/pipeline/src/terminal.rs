use std::path::PathBuf;

// TBDuse std::path::PathBuf;
use clap::{Args, Parser, Subcommand};

/// Cerebro: production stack server
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "cerebro-workflow", version)]
pub struct App {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Parse and process pipeline results
    Process(PipelineProcessArgs),
    
    /// Quality control table from processed results
    Quality(PipelineQualityArgs),
    /// Taxa table from processed results
    Taxa(PipelineTaxaArgs),
    /// Create a sample sheet from the input directory
    SampleSheet(PipelineSampleSheetArgs),

    #[clap(subcommand)]
    /// Utility tools for pipeline and stack
    Tools(ToolCommands),
}


#[derive(Debug, Subcommand)]
pub enum ToolCommands {
    /// Watch a file path with the provided configuration
    SubsetFasta(SubsetArgs),
    /// Send a message to the provided Slack channel
    Anonymize(AnonymizeArgs),
    /// Send a message to the provided Slack channel
    SplitFasta(SplitArgs),
    /// Watch a file path with the provided configuration
    ScanRemap(VirusScanRemapArgs),
    /// Send a message to the provided Slack channel
    UmiCheck(UmiCheckArgs),
    /// Send a message to the provided Slack channel
    UmiTrim(UmiTrimArgs),
    /// Send a message to the provided Slack channel
    UmiPrepCalib(UmiPrepareCalibArgs),
    /// Send a message to the provided Slack channel
    UmiDedupCalib(UmiDedupCalibArgs),
    /// Send a message to the provided Slack channel
    UmiDedupNaive(UmiDedupNaiveArgs),
}


#[derive(Debug, Args)]
pub struct SubsetArgs {
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
pub struct AnonymizeArgs {
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
pub struct SplitArgs {
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
===============================
PIPELINE PARSERS AND PROCESSORS
===============================
*/

#[derive(Debug, Args)]
pub struct PipelineProcessArgs {
    /// Pipeline sample results directories for processing
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Taxonomy directory containing 'names.dmp' and 'nodes.dmp' used for classification (NCBI)
    /// 
    /// Must be present for taxonomic processing, otherwise only the quality control module is proccessed.
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
pub struct PipelineTaxaArgs {
    /// Processed pipeline result samples (.json)
    #[clap(long, short = 'i', num_args(0..))]
    pub input: Vec<PathBuf>,
    /// Taxon filter configuration (.json)
    #[clap(long, short = 'f')]
    pub filter: Option<PathBuf>,
    /// Minimum total reads per million 
    #[clap(long, short = 'r', default_value="0")]
    pub min_rpm: f64,
    /// Minimum k-mer reads per million 
    #[clap(long, short = 'k', default_value="0")]
    pub min_rpm_kmer: f64,
    /// Minimum alignment reads per million 
    #[clap(long, short = 'a', default_value="0")]
    pub min_rpm_alignment: f64,
    /// Minimum alignment reads per million 
    #[clap(long, short = 'm', default_value="0")]
    pub min_rpm_remap: f64,
    /// Minimum number of assembled contigs
    #[clap(long, short = 'c', default_value="0")]
    pub min_contigs: u64,
    /// Minimum number of assembled bases
    #[clap(long, short = 'b', default_value="0")]
    pub min_bases: u64,
    /// Output file for taxa table
    #[clap(long, short = 'o')]
    pub output: PathBuf,
    /// Output delimited for taxa table
    #[clap(long, short = 's', default_value="\t")]
    pub sep: char,
    /// Extract sample identifier and tags for post-processing utilities
    #[clap(long, short = 'e',)]
    pub extract: bool,
    /// Include header in taxa table
    #[clap(long, short = 'H')]
    pub header: bool,

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
    pub symlinks: bool,
    /// Single reads for Oxford Nanopore Technology 
    #[clap(long)]
    pub ont: bool
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
