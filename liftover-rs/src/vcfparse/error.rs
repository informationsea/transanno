use nom::Offset;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum VCFParseError {
    #[error("I/O Error")]
    IoError(#[from] std::io::Error),
    #[error("Header parse error at line: {line}  column: {column}")]
    HeaderParseError { line: u32, column: usize },
    #[error("Record parse error at line: {line}  column: {column}")]
    RecordParseError { line: u32, column: usize },
    #[error("Allele frequency should be number at line: {line}")]
    FrequencyIsNotNumber { line: u32 },
    #[error("UTF-8 Error")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("multiple GT record")]
    InvalidGTRecord,
    #[error("parse float error")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("parse int error")]
    ParseIntError(#[from] std::num::ParseIntError),
    #[error("parse error")]
    ParseError,
}

impl From<nom::Err<(&[u8], nom::error::ErrorKind)>> for VCFParseError {
    fn from(_e: nom::Err<(&[u8], nom::error::ErrorKind)>) -> VCFParseError {
        VCFParseError::ParseError
    }
}

pub fn as_header_error(line: u32, data: &[u8], e: (&[u8], nom::error::ErrorKind)) -> VCFParseError {
    let column = data.offset(e.0);
    VCFParseError::HeaderParseError { line, column }.into()
}

pub fn as_record_error(line: u32, data: &[u8], e: (&[u8], nom::error::ErrorKind)) -> VCFParseError {
    let column = data.offset(e.0);
    VCFParseError::RecordParseError { line, column }.into()
}
