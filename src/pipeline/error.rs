/*
========================
Custom error definitions
========================
*/

use thiserror::Error;

#[derive(Error, Debug)]
pub enum WorkflowError {
    /// Indicates failure to parse a the fastp JSON
    #[error("failed to parse fastp output from JSON")]
    ParseFastp(#[source] serde_json::Error),
    /// Indicates failure to parse the Scrubby JSON
    #[error("failed to parse Scrubby output from JSON")]
    ParseScrubby(#[source] serde_json::Error),
    /// Represents all other cases of `std::io::Error`.
    #[error(transparent)]
    IOError(#[from] std::io::Error),
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
    TaxonomyNotParsed(#[source] taxonomy::Error)
}

#[derive(Error, Debug)]
pub enum WorkflowUtilityError {
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
    /// Represents all other cases of `csv::Error`.
    #[error("failed to process CSV")]
    CsvError(#[from] csv::Error),
    /// Represents all other cases of `std::io::Error`.
    #[error("failed to process input/output (IO)")]
    IoError(#[from] std::io::Error)
}
