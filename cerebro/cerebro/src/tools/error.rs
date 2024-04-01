use thiserror::Error;


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

    #[error("failed to hash password")]
    PasswordNotHashed,
    /// Indicates a failure to find the tax lineage for the provided taxid from the taxonomy
    #[error("could not find taxonomic name for taxid {1}")]
    TaxNameNotAvailable(#[source] taxonomy::Error, String)
}
