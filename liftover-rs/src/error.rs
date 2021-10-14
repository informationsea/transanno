use thiserror::Error;

#[derive(Debug, Error)]
pub enum LiftOverError {
    #[error("I/O Error")]
    IoError(#[from] std::io::Error),
    #[error("CSV Error")]
    CsvError(#[from] csv::Error),
    #[error("JSON Error")]
    JsonError(#[from] serde_json::Error),
    #[error("UTF-8 or IO Error")]
    Utf8IoError(#[from] std::str::Utf8Error),
    #[error("VCF parse error")]
    VCFParseError(#[from] crate::vcfparse::VCFParseError),
    #[error("Invalid number of header at line {0}")]
    InvalidNumberOfHeader(u32),
    #[error("Invalid header at line {0}. A header line should starts with \"chain\"")]
    NoChainHeaderFound(u32),
    #[error("Invalid strand at line {0}")]
    InvalidStrand(u32),
    #[error("Invalid chromosome length at line {0}")]
    InvalidChromosomeLength(u32),
    #[error("Unmatched chromosome length: {0} / chain file: {1} / FASTA: {2}")]
    UnmatchedChromosomeLength(String, u64, u64),
    #[error("Unmatched chromosome length: {0} / chain file: {1} / Reference FASTA: {2}")]
    UnmatchedReferenceChromosomeLength(String, u64, u64),
    #[error("Unmatched chromosome length: {0} / chain file: {1} / Query FASTA: {2}")]
    UnmatchedQueryChromosomeLength(String, u64, u64),
    #[error("Invalid number of columns at line {0}")]
    InvalidNumberOfColumns(u32),
    #[error("Reference sequence is different from expected sequence")]
    DifferentReference,
    #[error("Parse integer error")]
    ParseIntError(#[from] std::num::ParseIntError),
    #[error("Parse strand error")]
    ParseStrandError,
    #[error("stand for reference should be forward")]
    ReferenceStrandShouldForward,
    #[error(
        "length of chromosome {0} is not equal to length in chain file. Are you using correct reference?",
    )]
    ReferenceChromosomeLengthIsNotMatch(String),
    #[error(
         "length of chromosome {0} is not equal to length in chain file. Are you using correct query?",
    )]
    QueryChromosomeLengthIsNotMatch(String),
    #[error("unknown sequence error: {0}: {1}")]
    UnknownSequenceError(String, std::io::Error),
    #[error("Failed to parse gene annotation")]
    GeneParseError(#[from] crate::geneparse::GeneParseError),
    #[error("Chromosome {0} is not found in FASTA")]
    ChromosomeNotFound(String),
}
