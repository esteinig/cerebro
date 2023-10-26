use std::ffi::{OsString, OsStr};
use std::path::PathBuf;
use structopt::StructOpt;
use structopt::clap::AppSettings;
use thiserror::Error;


/// A collection of custom errors relating to the command line interface for this package.
#[derive(Error, Debug, PartialEq)]
pub enum CliError {
    /// Indicates that a string cannot be parsed into a [`CompressionFormat`](#compressionformat).
    #[error("{0} is not a valid output format")]
    CompressionFormat(String),
    /// Indicates that a string cannot be parsed into a [`CompressionLevel`](#compressionlevel).
    #[error("{0} is not a valid compression level (1-9)")]
    CompressionLevel(String),
    /// Indicates a bad combination of input and output files was passed.
    #[error("Bad combination of input and output files: {0}")]
    BadInputOutputCombination(String),
}

/// MGP-CEREBRO command-line application
#[derive(Debug, StructOpt)]
pub struct Cli {
    /// API Token File
    #[structopt(
        long = "--token-file"
    )]
    pub token_file: Option<PathBuf>,
    /// API Token ENV
    #[structopt(
        long = "--token-env",
        default_value = "CEREBRO_API_TOKEN"
    )]
    pub token_env: String,
    /// API URL
    #[structopt(
        long,
        default_value = "http://localhost:8080"
    )]
    pub url: String,
    /// DANGER: ignore SSL certificate verification!
    #[structopt(
        long
    )]
    pub danger_accept_invalid_certs: bool,
    #[structopt(subcommand)]
    pub commands: Commands
}


#[derive(Debug, StructOpt)]
pub enum Commands {
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    /// Loging user and return access token
    Login {
        /// Email for login
        #[structopt(
            short = "e",
            long,
        )]
        email: String,
        /// Password for login
        #[structopt(
            short = "p",
            long,
        )]
        password: Option<String>,
    },
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    /// Generate a report from a configuration
    Report {
        /// Output template file for rendering with `tectonic`
        #[structopt(
            short = "o",
            long,
        )]
        pdf: PathBuf,
        /// Output filled in template file for debugging
        #[structopt(
            short = "t",
            long,
        )]
        template: Option<PathBuf>,
        /// Base configuration file (TOML)
        #[structopt(
            short = "c",
            long,
        )]
        base_config: PathBuf,
        /// Full patient configuration file (TOML)
        #[structopt(
            short = "p",
            long,
        )]
        patient_config: Option<PathBuf>,
        /// Patient header configuration file (TOML)
        #[structopt(
            short = "H",
            long,
        )]
        patient_header_config: Option<PathBuf>,
        /// Patient header configuration file (TOML)
        #[structopt(
            short = "R",
            long,
        )]
        patient_result_config: Option<PathBuf>,
        /// Patient header configuration file (TOML)
        #[structopt(
            short = "s",
            long,
        )]
        sample_config: Option<Vec<PathBuf>>,
    },
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    /// Print a quick quality control summary to the console
    QcSummary {
        /// Paths to parsed workflow sample files (JSON)
        #[structopt(
            short = "i",
            long,
        )]
        input: Vec<PathBuf>,
        /// Path to QC table output CSV
        #[structopt(
            short = "o",
            long,
        )]
        output: Option<PathBuf>,
        /// Include header for output
        #[structopt(
            short = "H",
            long,
        )]
        header: bool
    },
    /// Run contaminant detection on a set of models
    Contam {
        /// Paths to Cerebro model files (JSON)
        #[structopt(
            short = "i",
            long,
        )]
        input: Vec<PathBuf>,
        /// Include samples with these tags 
        #[structopt(
            short = "y",
            long,
        )]
        include_tags: Vec<String>,
        /// Exclude samples with these tags 
        #[structopt(
            short = "n",
            long,
        )]
        exclude_tags: Vec<String>,
        /// Minimum ERCC constructs to pass filter
        #[structopt(
            short = "e",
            long,
            default_value = "60"
        )]
        min_ercc: usize, 
        /// Maximum low completexity percent of total reads
        #[structopt(
            short = "c",
            long,
            default_value = "20"
        )]
        max_low_complexity: f32, 
        /// Minimum percent of samples in filtered set
        /// that a taxon must occurr on to be considered
        /// as potential contaminant
        #[structopt(
            short = "s",
            long,
            default_value = "0"
        )]
        samples: f32,
        /// Minimum percent of samples in filtered set
        /// that a taxon must occurr on to be considered
        /// as potential contaminant
        #[structopt(
            short = "d",
            long
        )]
        domains: Vec<String>,
        /// Path to summary table output (CSV)
        #[structopt(
            short = "o",
            long,
        )]
        output: Option<PathBuf>,
        /// ERCC input mass in picogram
        #[structopt(
            short = "m",
            long,
            default_value = "25"
        )]
        ercc_mass: u32,
        /// ERCC input mass in picogram
        #[structopt(
            short = "t",
            long,
        )]
        taxid_mass: Option<Vec<String>>
    },
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    /// Parse a directory of paired-end read files to create a sample sheet
    CreateSampleSheet {
        /// Directory containing the paired-end input read files
        #[structopt(
            short = "i",
            long,
            parse(try_from_os_str = check_file_exists),
            required = true,
            multiple = true,
        )]
        input: Vec<PathBuf>,
        /// Output sample sheet (CSV)
        #[structopt(
            short = "o",
            long,
            required = true
        )]
        output: PathBuf,
        /// Sample glob - pattern matching to find paired samples
        /// 
        /// The glob string to specify paired sample extensions should be in format: {forward_extension,reverse_extension}
        /// where the wildcard specifies the sample identifier, for example: *_{R1,R2}.fastq.gz will match files
        /// "Sample1_R1.fastq.gz" and "Sample1_R2.fastq.gz" to sample identifier "Sample1"
        #[structopt(
            short = "g",
            long,
            default_value = "*_{R1_001,R2_001}.fastq.gz"
        )]
        glob: String,
        /// Run identifier - if not provided uses input directory name
        /// 
        /// If you want to fill in custom run identifiers for each sample, 
        /// you can provide an empty string ("") and edit the sample sheet 
        /// manually.
        #[structopt(
            short = "r",
            long,
        )]
        run_id: Option<String>,
        /// Run date - if not provided uses current date (YYYYMMDD)
        /// 
        /// If you want to fill in custom run dates for each sample, 
        /// you can provide an empty string ("") and edit the sample sheet 
        /// manually.
        #[structopt(
            short = "d",
            long,
        )]
        run_date: Option<String>,
        /// Sample group - if not provided group designation is an empty string
        /// 
        /// Sample groups can be specified manually for larger runs containing
        /// sampels from multiple experimental groups - these are later available
        /// in the front-end application.
        #[structopt(
            short = "s",
            long,
        )]
        sample_group: Option<String>,
        /// Allow symlink target reading for glob file walking
        #[structopt(
            short = "l",
            long,
        )]
        symlinks: bool
        
    },
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    /// Parse workflow data for a single sample
    ParseSample {
        /// Results directory for a sample from workflow output.
        #[structopt(
            short = "i",
            long,
            parse(try_from_os_str = check_file_exists),
            required = true
        )]
        input: PathBuf,
        /// Sample identifier to use for output, otherwise input path stem.
        #[structopt(
            short = "s",
            long
        )]
        sample_id: Option<String>,
        /// Directory containing `names.dmp` and `nodes.dmp`.
        #[structopt(
            short = "n",
            long,
            required = true
        )]
        taxonomy: PathBuf,
        /// Output filename for sample data (JSON)
        #[structopt(
            short = "o",
            long
        )]
        output: Option<PathBuf>,
        /// Output all aggregated taxon evidence
        /// 
        /// By default we aggregate evidence at species level, that
        /// is evidence from below species level is pulled up the
        /// lineage to species level, and any taxa above species
        /// level are discarded.
        #[structopt(
            short = "a",
            long
        )]
        all_taxa: bool,
    },
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    /// Read a sample from JSON and insert into the DB
    InsertSample {
        /// JSON file of parsed sample.
        #[structopt(
            short = "i",
            long,
            parse(try_from_os_str = check_file_exists),
            required = true
        )]
        input: PathBuf,
        /// Workflow sample sheet to match run identifier to sample.
        #[structopt(
            short = "s",
            long,
            parse(try_from_os_str = check_file_exists),
            required = true
        )]
        sample_sheet: PathBuf,
        /// Workflow configuration file for database model.
        #[structopt(
            short = "w",
            long,
            parse(try_from_os_str = check_file_exists),
            required = true
        )]
        workflow_config: PathBuf,
        /// Team name of project collection to insert model into.
        #[structopt(
            short = "t",
            long
        )]
        team_name: String,
        /// Project name to insert model into.
        #[structopt(
            short = "p",
            long
        )]
        project_name: String,
        /// Database name to insert model into (explicit)
        #[structopt(
            short = "d",
            long
        )]
        db_name: Option<String>,
        /// Output model as JSON.
        #[structopt(
            short = "o",
            long,
        )]
        output: Option<PathBuf>,
    },
    /// Subset a database sequence file
    Subset {
        /// Database sequence file to subset (FASTA)
        #[structopt(
            short, long, parse(try_from_os_str = check_file_exists), required = true
        )]
        fasta: PathBuf,
        /// Mash screen output for subsetting (can be winner-take-all)
        #[structopt(
            short, long, parse(try_from_os_str = check_file_exists)
        )]
        mash: Option<PathBuf>,
        /// Minimum query identity to include sequence in subset [Mash]
        #[structopt(short = "i", long, default_value = "0")]
        min_identity: f64,
        /// Minimum number of shared hashes to include sequence in subset [Mash]
        #[structopt(short = "s", long, default_value = "0")]
        min_shared_hashes: u32,
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
        #[structopt(short = "g", long)]
        group_index: Option<usize>,
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
        #[structopt(short = "", long)]
        group_by: Option<String>,
        /// Field delimiter in identifier of mash query sequence 
        #[structopt(short = "d", long, default_value = "|")]
        group_sep: String,
        /// Output filepath with sequences removed
        #[structopt(short, long, parse(from_os_str), required = true)]
        output: PathBuf,
        /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
        ///
        /// Default is to attempt to infer the output compression format automatically from the filename
        /// extension (gz|bz|bz2|lzma). This option is used to override that.
        #[structopt(
            short = "O",
            long,
            value_name = "u|b|g|l",
            parse(try_from_str = parse_compression_format),
            possible_values = &["u", "b", "g", "l"],
            case_insensitive = true,
            hide_possible_values = true
        )]
        output_format: Option<niffler::compression::Format>,
        /// Compression level to use if compressing output
        #[structopt(
            short = "l",
            long,
            parse(try_from_str = parse_level),
            default_value = "6",
            value_name = "1-9"
        )]
        compression_level: niffler::Level,
    },
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    Completeness {
        /// Input filepath(s) (fa, gz, bz) for FASTA.
        ///
        /// FASTA file to compute completeness (100-percent of missing characters) for each sequence therein contained.
        /// Counts all "N" characters and divides byt the total length of the sequence and outputs to stdout the 
        /// following tab-delimited fields: seq id, n count, seq length, completeness
        #[structopt(
            short,
            long,
            required = true
        )]
        input: PathBuf,
        /// Output filepath(s) for (fa, gz, bz) for FASTA.
        #[structopt(
            short, 
            long, 
            parse(from_os_str), 
            required = true
        )]
        output: PathBuf,
        /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
        ///
        /// Default is to attempt to infer the output compression format automatically from the filename
        /// extension (gz|bz|bz2|lzma). This option is used to override that.
        #[structopt(
            short = "0",
            long,
            value_name = "u|b|g|l",
            parse(try_from_str = parse_compression_format),
            possible_values = &["u", "b", "g", "l"],
            case_insensitive=true,
            hide_possible_values = true
        )]
        output_format: Option<niffler::compression::Format>,
        /// Compression level to use if compressing output
        #[structopt(
            short = "l",
            long,
            parse(try_from_str = parse_level),
            default_value = "6",
            value_name = "1-9"
        )]
        compression_level: niffler::Level,
        /// Server port
        ///
        #[structopt(
            short = "m",
            long,
            required = false,
            default_value = "95.0"
        )]
        min_completeness: f64,
    },
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    Annotate {
        /// Input filepath(s) (fa, gz, bz) for database to annotate.
        ///
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists),
            required = true
        )]
        input: PathBuf,
        /// Output filepath(s) (fa, gz, bz) for annotated database.
        ///
        #[structopt(
            short,
            long,
            required = true
        )]
        output: PathBuf,
        /// Input filepath(s)of uncompressed accession to taxid file for annotation.
        ///
        #[structopt(
            short,
            long,
            multiple = true,
            required = true
        )]
        acc2taxid: Vec<PathBuf>,
        /// Directory containing `names.dmp` and `nodes.dmp`
        #[structopt(
            short,
            long
        )]
        ncbi_tax: Option<PathBuf>,
        /// Directory containing `names.dmp` and `nodes.dmp`
        #[structopt(
            short,
            long
        )]
        exclude: Vec<String>,
    },
    /// Get the taxonomy metadata for `mgp-cerebro`
    GetDbMetadata {
        /// Workign directory for downloads and database building
        ///
        #[structopt(
            short = "w",
            long
        )]
        workdir: Option<PathBuf>,
        /// Force use working directory if it already exists
        ///
        #[structopt(
            short = "f",
            long
        )]
        force: bool
    },
    /// Get the genome sequences for database auto-builds
    GetDbData {
        /// Workign directory for downloads and database building
        ///
        #[structopt(
            short = "w",
            long
        )]
        workdir: Option<PathBuf>,
        /// Force use working directory if it already exists
        ///
        #[structopt(
            short = "f",
            long
        )]
        force: bool
    },
    /// Run the server for the API
    RunServer {
        /// Server host
        ///
        #[structopt(
            short = "s",
            long,
            required = true,
            default_value = "127.0.0.1"
        )]
        host: String,
        /// Server port
        ///
        #[structopt(
            short = "p",
            long,
            required = true,
            default_value = "8080"
        )]
        port: u16,
        /// Threads for server
        ///
        #[structopt(
            short = "t",
            long,
            default_value = "4"
        )]
        threads: usize,
        /// Configuration file
        ///
        #[structopt(
            short = "c",
            long,
        )]
        toml: Option<PathBuf>,
        /// Configuration file env var
        ///
        #[structopt(
            short = "e",
            long,
        )]
        env: Option<String>,
    },
    #[structopt(global_settings = &[AppSettings::ColoredHelp])]
    PingApi { },
    /// Anonymize read files
    Anonymize {
        /// Input filepath(s) (fa, fq, gz, bz).
        ///
        /// For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two
        /// files consecutively `-i r1.fq r2.fq`. NOTE: Read identifiers for paired-end Illumina reads
        /// are assumed to be the same in forward and reverse read files (modern format) without trailing
        /// read orientations e.g. `/1` and `/2`. If you are using legacy identifiers, reads in the anonymized
        /// output may be unpaired.
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists),
            multiple = true,
            required = true
        )]
        input: Vec<PathBuf>,
        /// Output filepath(s) with anonymized read identifiers.
        ///
        /// For paired Illumina you may either pass this flag twice `-o r1.fq -o r2.fq` or give two
        /// files consecutively `-o r1.fq r2.fq`. NOTE: The order of the pairs is assumed to be the
        /// same as that given for --input.
        #[structopt(short, long, parse(from_os_str), multiple = true, required = true)]
        output: Vec<PathBuf>,
        /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
        ///
        /// Default is to attempt to infer the output compression format automatically from the filename
        /// extension (gz|bz|bz2|lzma). This option is used to override that.
        #[structopt(
            short = "O",
            long,
            value_name = "u|b|g|l",
            parse(try_from_str = parse_compression_format),
            possible_values = &["u", "b", "g", "l"],
            case_insensitive=true,
            hide_possible_values = true
        )]
        output_format: Option<niffler::compression::Format>,
        /// Compression level to use if compressing output
        #[structopt(
            short = "l",
            long,
            parse(try_from_str = parse_level),
            default_value = "6",
            value_name = "1-9"
        )]
        compression_level: niffler::Level,
        /// Add a fake Illumina header instead of UUID (paired end data only)
        #[structopt(
        short = "-f",
        long
        )]
        fake_illumina_header: bool 
    },
    Split {
        /// Input filepath(s) (fa, gz, bz) for FASTA file to split.
        ///
        /// FASTA file to split. Each sequence in the FASTA file is 
        /// prefixed with the index of the sequence in the input file
        /// and the input filename (e.g. {idx}-{fname}.fasta).
        /// 
        /// Multiple FASTA files can be provided to split. You may either 
        /// pass this flag twice `-i f1.fasta -i f2.fasta` or give two files
        /// consecutively `-i f1.fasta f2.fasta`
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists),
            multiple = true,
            required = true
        )]
        input: Vec<PathBuf>,
        /// Output directory that contains the split fasta files
        /// 
        /// If multiple input files are given, split output files are
        /// indexed anew for each file with the seqeunce index and the 
        /// file name. Output compression is either inferred from the filename
        /// or specified can be specified --output-format (gz|bz|bz2|lzma)
        #[structopt(
            short,
            long,
            required = true
        )]
        outdir: PathBuf,
        /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
        ///
        /// Default is to attempt to infer the output compression format automatically from the filename
        /// extension (gz|bz|bz2|lzma). This option is used to override that.
        #[structopt(
            short = "O",
            long,
            value_name = "u|b|g|l",
            parse(try_from_str = parse_compression_format),
            possible_values = &["u", "b", "g", "l"],
            case_insensitive=true,
            hide_possible_values = true
        )]
        output_format: Option<niffler::compression::Format>,
        /// Compression level to use if compressing output
        #[structopt(
            short = "l",
            long,
            parse(try_from_str = parse_level),
            default_value = "6",
            value_name = "1-9"
        )]
        compression_level: niffler::Level
    },
    HashPassword {
        /// Password string to hash
        #[structopt(
            short,
            long,
            required = true
        )]
        password: String,
    }
}


// Functions may be heavily adapted from Rasusa, due to the excellent error annotation style
impl Cli {
    /// Checks there is a valid and equal number of `--input` and `--output` arguments given.
    ///
    /// # Errors
    /// 
    /// A [`CliError::BadInputOutputCombination`](#clierror) is returned for the following:
    /// - Either `--input` or `--output` are passed more than twice
    /// - An unequal number of `--input` and `--output` are passed
    pub fn validate_input_output_combination(&self) -> Result<(), CliError> {
        match &self.commands {
            Commands::Anonymize { input, output, .. } => {
                let out_len = output.len();
                let in_len = input.len();
                if in_len > 2 {
                    let msg = String::from("Got more than 2 files for input.");
                    return Err(CliError::BadInputOutputCombination(msg));
                }
                if out_len > 2 {
                    let msg = String::from("Got more than 2 files for output.");
                    return Err(CliError::BadInputOutputCombination(msg));
                }
                if in_len != out_len {
                    let msg = format!("Got {} --input but {} --output", in_len, out_len);
                    return Err(CliError::BadInputOutputCombination(msg));
                }
            },
            _ => {},
        };
        Ok(())
    }
}

/// A utility function to validate compression format is in allowed values
fn parse_compression_format(s: &str) -> Result<niffler::compression::Format, CliError> {
    match s {
        "b" | "B" => Ok(niffler::Format::Bzip),
        "g" | "G" => Ok(niffler::Format::Gzip),
        "l" | "L" => Ok(niffler::Format::Lzma),
        "u" | "U" => Ok(niffler::Format::No),
        _ => Err(CliError::CompressionFormat(s.to_string())),
    }
}

/// A utility function to validate compression level is in allowed range
#[allow(clippy::redundant_clone)]
fn parse_level(s: &str) -> Result<niffler::Level, CliError> {
    let lvl = match s.parse::<u8>() {
        Ok(1) => niffler::Level::One,
        Ok(2) => niffler::Level::Two,
        Ok(3) => niffler::Level::Three,
        Ok(4) => niffler::Level::Four,
        Ok(5) => niffler::Level::Five,
        Ok(6) => niffler::Level::Six,
        Ok(7) => niffler::Level::Seven,
        Ok(8) => niffler::Level::Eight,
        Ok(9) => niffler::Level::Nine,
        _ => return Err(CliError::CompressionLevel(s.to_string())),
    };
    Ok(lvl)
}



/// A utility function to validate whether an input files exist
fn check_file_exists(file: &OsStr) -> Result<PathBuf, OsString> {
    let path = PathBuf::from(file);
    let path_msg = format!("{:?} does not exist", path);
    if path.exists() {
        let abs_path = path.canonicalize().map_err(|_| OsString::from(path_msg))?;
        Ok(abs_path)
    } else {
        Err(OsString::from(path_msg))
    }
}
