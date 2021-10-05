use nom::error::ErrorKind;
use nom::Err;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum GeneParseError {
    #[error("I/O Error")]
    IoError(#[from] std::io::Error),
    #[error("Parse error")]
    ParseError,
    #[error("Parse int error")]
    ParseIntError(#[from] std::num::ParseIntError),
    #[error("Parse error")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("Parse error at line: {}", _0)]
    ParseErrorAtLine(u64, Box<GeneParseError>),
    #[error("Gene grouping error")]
    GroupingError,
    #[error("# of refGene columns should be 16")]
    RefGeneInvalidColumnNumber,
    #[error("Invalid strand")]
    StrandParseError,
    #[error("CDS stat parse error")]
    CDSStatParseError,
}

impl From<Err<(&str, ErrorKind)>> for GeneParseError {
    fn from(e: Err<(&str, ErrorKind)>) -> GeneParseError {
        println!("nom error: {:?}", e);
        GeneParseError::ParseError
    }
}
