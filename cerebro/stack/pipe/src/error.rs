/*
========================
Custom error definitions
========================
*/

use thiserror::Error;

#[derive(Error, Debug)]
pub enum WorkflowError {

    #[error(transparent)]
    CsvError(#[from] csv::Error),
    #[error(transparent)]
    SerdeJsonError(#[from] serde_json::Error),
    #[error(transparent)]
    ScrubbyError(#[from] scrubby::error::ScrubbyError),
    #[error(transparent)]
    VircovError(#[from] vircov::error::VircovError),
    /// Represents all other cases of `reqwest::Error`.
    #[error(transparent)]
    ReqwestError(#[from] reqwest::Error),
    #[error(transparent)]
    IoError(#[from] std::io::Error),
    /// Represents errors from building a Rayon thread pool.
    #[error(transparent)]
    RayonThreadPoolError(#[from] rayon::ThreadPoolBuildError),
    /// Represents all other cases of `niffler::Error`.
    #[error(transparent)]
    NifflerError(#[from] niffler::Error),

    /// Failed to make the download request
    #[error("failed to execute request: {0}")]
    DownloadFailedRequest(reqwest::StatusCode),
    /// Failed to configure the downloader through the builder pattern due to missing field
    #[error("failed to configure the output directory field for the downloader")]
    DownloaderMissingOutdir,
    
    #[error("Failed to convert OsString to String")]
    FileNameConversionError,

    #[error("Failed to detect required output file")]
    PipelineOutputNotFound,

    /// Represents an error when failing to extract a sequence record header.
    #[error("failed to extract sequence record header")]
    NeedletailHeader(#[source] std::str::Utf8Error),
    /// Represents an error when failing to extract a valid header of a read.
    #[error("failed to extract a valid header of read")]
    NeedletailFastqHeader,
    
    /// Indicates failure to parse a the nanoq JSON
    #[error("failed to parse nanoq output from JSON")]
    ParseNanoq(#[source] serde_json::Error),
    /// Indicates failure to parse the Scrubby JSON
    #[error("failed to parse Scrubby output from JSON")]
    ParseScrubby(#[source] serde_json::Error),
    /// Indicates failure to find the required nanoq scan output
    #[error("failed to find nanoq scan output")]
    NanoqScan,
    /// Indicates failure to find a taxid in the scan-remap record and subsequent failure to aggregate the record
    #[error("failed to find taxonomic identifier in Vircov Scan-Remap Record")]
    VircovScanRemapTaxidMissing,
    /// Represents a failure to convert the field from string to float field
    #[error("failed to convert a float field from `BLAST LCA`")]
    BlastLcaFloatFieldConversion(#[source] std::num::ParseFloatError),
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the read field in the report from `Kraken2Uniq`")]
    KrakenReportReadFieldConversion,
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the direct read field in the report from `Kraken2Uniq`")]
    KrakenReportDirectReadFieldConversion,
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the kmer field in the report from `Kraken2Uniq`")]
    KrakenReportKmerFieldConversion,
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the unique kmer field in the report from `Kraken2Uniq`")]
    KrakenReportKmerUniqueFieldConversion,
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the regions field in the output from `Vircov`")]
    VircovRegionFieldConversion,
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the reads field in the output from `Vircov`")]
    VircovReadFieldConversion,
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the alignments field in the output from `Vircov`")]
    VircovAlignmentFieldConversion,
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the basepair field in the output from `Vircov`")]
    VircovBasepairFieldConversion,
    /// Represents a failure to convert the field from string to numeric field in the report file
    #[error("failed to convert the reference length field in the output from `Vircov`")]
    VircovLengthFieldConversion,
    /// Represents a failure to convert the field from string to float field in the report file
    #[error("failed to convert the coverage field in the output from `Vircov`")]
    VircovCoverageFieldConversion,
    /// Represents a failure to convert a PathBuf converted OsStr into a String
    #[error("incorrect format of the reference path, are there non-UTF8 characters?")]
    InvalidReferencePath,
    /// Indicates a failure to obtain an absolute path
    #[error("reference name could not be obtained from {0}")]
    ReferenceNameExtraction(String),
    /// Indicates a failure to find the index in a result file name
    #[error("reference name could not be obtained from {0} at index {1}")]
    ResultFileIndexNotFound(String, String),
    /// Indicates a failure to find a specific tax rank for the provided taxid in the taxonomy
    #[error("could not find taxonomic rank `(1)` for taxid `{2}`")]
    TaxRankNotFound(#[source] taxonomy::Error, String, String),
    /// Indicates a failure to rcover the tax rank for the provided taxid from the taxonomy
    #[error("could not find taxonomic rank for taxid {1}")]
    TaxRankNotAvailable(#[source] taxonomy::Error, String),
    /// Indicates a failure to find the tax lineage for the provided taxid from the taxonomy
    #[error("could not find taxonomic lineage for taxid {1}")]
    TaxLineageNotAvailable(#[source] taxonomy::Error, String),
    /// Indicates a failure to find the tax lineage for the provided taxid from the taxonomy
    #[error("could not find taxonomic name for taxid {1}")]
    TaxNameNotAvailable(#[source] taxonomy::Error, String),
    /// Indicates a failure to call a taxon aggregation functions again
    #[error("re-aggregation of taxa is not allowed for {0}")]
    TaxAggregate(String),
    /// Indicates a failure to build a glob from multiple patterns for file matching
    #[error("failed to build glob walker for file parsing")]
    GlobWalkBuilder,
    /// Indicates a failure to parse a traxid from an alignment output sequence description
    #[error("failed to parse the taxid field from a Vircov record header (DB: {0})")]
    VircovTaxidFieldMissing(String),
    /// Indicates failure with JSON serialization
    #[error("failed to serialize JSON")]
    JsonSerialization(#[source] serde_json::Error),
    /// Indicates a failure to obtain a row from the sample summary table of taxa
    #[error("failed to obtain a row from the sample summary of taxa")]
    SampleSummaryTableRows,
    // Indicates a failure to obtain a row from the sample summary quality control 
    #[error("failed to obtain a row from the sample summary of quality control data")]
    SampleSummaryQualityRows,
    /// Indicates a failure to obtain a record from the vector of BlastRecords
    #[error("failed to obtain a single record for lca record length extraction")]
    BlastLcaRecordExtraction,
    /// Indicates a failure to obtain a coverage value from a split contig name vector of a BlastRecord
    #[error("failed to obtain the coverage value from a vector of split name strings")]
    BlastLcaCovExtraction,
    /// Indicates a failure to obtain a taxids for LCA lineage computations
    #[error("failed to obtain some taxonomic identifiers for LCA computation")]
    BlastLcaTaxidsNotFound,
    /// Indicates a failure to obtain the record with the lowest evalue from a vector of BlastRecords 
    #[error("failed to obtain the lowest evalue record")]
    BlastLcaRecordExtractionEvalue,
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an integer field `query length` from `BLAST`")]
    BlastQueryLengthFieldConversion(#[source] std::num::ParseIntError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an integer field `query start` from `BLAST`")]
    BlastQueryStartFieldConversion(#[source] std::num::ParseIntError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an integer field `query end` from `BLAST`")]
    BlastQueryEndFieldConversion(#[source] std::num::ParseIntError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an integer field `subject length` from `BLAST`")]
    BlastSubjectLengthFieldConversion(#[source] std::num::ParseIntError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an integer field `subject start` from `BLAST`")]
    BlastSubjectStartFieldConversion(#[source] std::num::ParseIntError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an integer field `subject end` from `BLAST`")]
    BlastSubjectEndFieldConversion(#[source] std::num::ParseIntError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an integer field `length` from `BLAST`")]
    BlastLengthFieldConversion(#[source] std::num::ParseIntError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an integer field `identity number` from `BLAST`")]
    BlastIdentityNumberFieldConversion(#[source] std::num::ParseIntError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an float field `identity percent` from `BLAST`")]
    BlastIdentityPercentFieldConversion(#[source] std::num::ParseFloatError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an float field `evalue` from `BLAST`")]
    BlastEvalueFieldConversion(#[source] std::num::ParseFloatError),
    /// Represents a failure to convert the field from string to numeric field
    #[error("failed to convert an float field `bitscore` from `BLAST`")]
    BlastBitscoreFieldConversion(#[source] std::num::ParseFloatError),
    /// Represents a failure to summarize quality control data
    #[error("failed to summarize quality control data as quality control module is incomplete")]
    QualityControlNotConducted,
    /// Represents a failure to parse the reference taxonomy
    #[error("failed to parse the reference taxonomy")]
    TaxonomyNotParsed(#[source] taxonomy::Error),
    /// Represents a failure to erialize a row
    #[error("failed to serialize quality control table row")]
    SerializeQualityTableRow(#[source] csv::Error),
    /// Represents a failure to serialize a row
    #[error("failed to serialize taxa table row")]
    SerializeTaxaTableRow(#[source] csv::Error),
    /// Represents a failure to serialize a row
    #[error("failed to create the taxon output file")]
    CreateTaxaOutputFile(#[source] csv::Error),
    /// Indicates failure to parse scannign stage alignment summary
    #[error("failed get parse vircov records from scanning stage")]
    VircovScanParseFailed,
    /// Indicates failure to parse remap stage alignment summary
    #[error("failed get parse vircov records from remap stage")]
    VircovRemapParseFailed,
    /// Indicates failure to find a matching scan record for a remap record - 
    /// remap records should always subsets of scan records
    #[error("failed get a matching remap record for the scan record of reference: {0}")]
    VircovRemapMatchingScanRecordNotFound(String),
    /// Indicates failure to assign a consensus sequence to the alignment result,
    /// likely due to not accounting for segmentation (and aligning + building
    /// the consensus sequnces for each segment in the pipeline)
    #[error("multiple consensus sequences found for reference: {0}")]
    VircovConsensusSequenceMatchNotIdentifiable(String),
    /// Indicates failure to assign a remap result to a single scan result,
    /// likely due to not accounting for segmentation (and aligning + building
    /// the consensus sequnces for each segment in the pipeline)
    #[error("multiple remap references found for reference: {0}")]
    VircovRemapMatchNotIdentifiable(String),
    /// Indicates failure to assign a coverage result to a single scan result,
    /// likely due to not accounting for segmentation (and aligning + building
    /// the consensus sequnces for each segment in the pipeline)
    #[error("multiple coverage references found for reference: {0}")]
    VircovCoverageMatchNotIdentifiable(String),
    /// Indicates failure to deserialize a CSV record into struct
    #[error("failed deserialize a csv record")]
    CsvRecordNotParsed(#[source] csv::Error),
    /// Indicates failure to parse sequence record identifier from Needletail
    #[error("failed to parse sequence record identifier")]
    RecordIdentifierNotParsed,
    /// Indicates failure to parse a record name from Needletail
    #[error("failed to parse sequence record header")]
    RecordHeaderUtf8NotParsed(#[from] std::str::Utf8Error),
    /// Indicates failure to parse file with Needletail
    #[error("failed read file input")]
    FastxParser(#[from] needletail::errors::ParseError),
    /// Indicates failure to specify a matchign number of input and 
    /// output file paths used in methods with parallel file processing 
    #[error("number of input and output files does not match")]
    InputOutputFileMismatch,
    /// Indicates failure to parse the last field of the sequence identifier which shoudl contain the UMI
    #[error("failed to parse the last sequence identifier field with the given separator: {0}")]
    UmiFieldNotFound(String),
    /// Indicates failure to parse an integer in a record
    #[error("failed to parse a valid integer into record")]
    RecordIntError(#[from] std::num::ParseIntError),
    /// Indicates failure to parse a float in a record
    #[error("failed to parse a valid float into record")]
    RecordFloatError(#[from] std::num::ParseFloatError),
    /// Indicates failure to get quality values from fastq
    #[error("failed get quality scores from fastq")]
    QualityScoresNotPresent,
    /// Indicates failure to find a matching UMI record from the read identifer
    #[error("failed get unique molecular identifier record matching the read identifier: {0}")]
    UmiNotFound(String),
    /// Indicates failure to deserialize cluste record from Calib
    #[error("failed to deserialize cluster record")]
    CalibClusterRecordNotValid(#[source] csv::Error),
    /// Indicates failure to deduplicate a cluster from Calib
    #[error("failed to select a read identifier from a cluster based on best average paired quality")]
    CalibClusterPickNotAvailable,
    /// Indicates failure to parse a group by index in a record identifier
    #[error("failed to parse group index [{0}] from record identifier [{1}]")]
    GroupIndexError(String, String),
    /// Indicates failure to select an optimal record from a group
    #[error("failed to select a record for group [{0}] ")]
    GroupSelectionError(String),
    /// Indicates failure to extract file name
    #[error("failed to extract file name")]
    FileName,
    /// Indicates failure with JSON deserialization.
    #[error("failed to deserialize JSON")]
    JsonDeserialization(#[source] serde_json::Error),

    /// Indicates failure to compile the sample id regex
    #[error("failed to compile sample identifier regex")]
    SampleIdRegex(#[source] fancy_regex::Error),
    /// Indicates failure to obtain the sample id regex capture
    #[error("failed to get sample identifier regex capture")]
    SampleIdRegexCapture(#[source] fancy_regex::Error),
    /// Indicates failure to obtain the sample id regex capture match
    #[error("failed to get sample identifier regex capture match")]
    SampleIdRegexCaptureMatch,
    /// Indicates failure to compile the sample tag regex
    #[error("failed to compile sample tag regex")]
    SampleTagRegex(#[source] fancy_regex::Error),
    /// Indicates failure to obtain the sample tag regex capture
    #[error("failed to get sample tag regex capture")]
    SampleTagRegexCapture(#[source] fancy_regex::Error),
    /// Indicates failure to obtain the sample tag regex capture match
    #[error("failed to get sample tag regex capture match")]
    SampleTagRegexCaptureMatch,

    /// Represents a failure to extract the sample identifier from a glob matched pattern
    #[error("failed to extract pattern match for sample identifier from: {0}")]
    GlobMatchSampleIdentifier(String),
    /// Represents a failure to 
    #[error("failed to extract an entry from the globbed walk through directory: {0}")]
    GlobWalk(String),
    /// Represents a failure to 
    #[error("failed to create a glob matcher for pattern: {0}")]
    GlobCreate(String),
    /// Represents a failure to 
    #[error("failed to find paired files for sample: {0}")]
    GlobPairedFiles(String),
    /// Represents a failure to 
    #[error("failed to create a sample sheet entry, number of files is incorrect: {0}")]
    SampleSheetEntryFiles(String),
    /// Represents a failure to 
    #[error("sample identifiers in sample sheet are not unique")]
    SampleSheetIdsNotUnique,
    /// Represents a failure to 
    #[error("failed to create a sample sheet entry, number of files is incorrect: {0}")]
    SampleSheetEntryRunDirectoryPath(String),
    /// Represents a failure to initiate writer
    #[error("failed to create a sample sheet writer to file: {0}")]
    SampleSheetCsvWriter(String),
    /// Represents a failure to initiate reader
    #[error("failed to create a sample sheet reader from file: {0}")]
    SampleSheetCsvReader(String),

    /// Represents a failure to reover Sylph lineage data
    #[error("failed to recover tax info from Sylph output lineage: {0}")]
    SylphTaxInfoRecoveryFailure(String),
}

#[derive(Error, Debug)]
pub enum WorkflowUtilityError {
   
    /// Represents all other cases of `csv::Error`.
    #[error("failed to process CSV")]
    CsvError(#[from] csv::Error),
    /// Represents all other cases of `std::io::Error`.
    #[error("failed to process input/output (IO)")]
    IoError(#[from] std::io::Error)
}
