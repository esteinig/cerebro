use thiserror::Error;

use crate::pipeline::error::WorkflowError;

/*
========================
Custom error definitions
========================
*/
#[derive(Error, Debug)]
pub enum ToolError {

    /// Indicates failure input/output file
    #[error("failed to read file")]
    FileIO(#[from] std::io::Error),
    /// Indicates failure to specify a matchign number of input and 
    /// output file paths used in methods with parallel file processing 
    #[error("number of input and output files does not match")]
    InputOutputFileMismatch,
    /// Indicates failure to parse file with Needletail
    #[error("failed read file input")]
    FastxParser(#[from] needletail::errors::ParseError),
    /// Indicates failure to parse a record name from Needletail
    #[error("failed to parse sequence record header")]
    RecordHeaderUtf8NotParsed(#[from] std::str::Utf8Error),
    /// Indicates failure to parse sequence record identifier from Needletail
    #[error("failed to parse sequence record identifier")]
    RecordIdentifierNotParsed,
    /// Indicates failure to parse the last field of the sequence identifier which shoudl contain the UMI
    #[error("failed to parse the last sequence identifier field with the given separator: {0}")]
    UmiFieldNotFound(String),
    /// Indicates failure to obtain compression writer with Niffler
    #[error("failed to get compression writer")]
    CompressionWriter(#[from] niffler::Error),
    /// Indicates failure to parse an integer in a record
    #[error("failed to parse a valid integer into record")]
    RecordIntError(#[from] std::num::ParseIntError),
    /// Indicates failure to parse a float in a record
    #[error("failed to parse a valid float into record")]
    RecordFloatError(#[from] std::num::ParseFloatError),

    /// Indicates failure to get quality values from fastq
    #[error("failed get quality scores from fastq")]
    QualityScoresNotPresent,

    /// Indicates failure to process a Vircov Report
    #[error("failed process vircov report")]
    VircovError(#[source] WorkflowError),


    /// Indicates failure to deserialize a CSV record into struct
    #[error("failed deserialize a csv record")]
    CsvRecordNotParsed(#[source] csv::Error),

    /// Failed to extract file stem of input file
    #[error("failed get extract file stem")]
    FileStemNotExtracted(#[source] WorkflowError),

    /// Indicates failure to find a matching UMI record from the read identifer
    #[error("failed get unique molecular identifier record matching the read identifier: {0}")]
    UmiNotFound(String),


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


    #[error("failed to hash password")]
    PasswordNotHashed,


    /// Indicates failure to extract file name
    #[error("failed to extract file name")]
    FileName,


    /// Indicates a failure to find the tax lineage for the provided taxid from the taxonomy
    #[error("could not find taxonomic name for taxid {1}")]
    TaxNameNotAvailable(#[source] taxonomy::Error, String)
}
