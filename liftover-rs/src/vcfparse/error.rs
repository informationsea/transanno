use failure::*;
use nom::Offset;
use std::fmt::{self, Display};

#[derive(Debug)]
pub struct VCFParseError {
    inner: failure::Context<VCFParseErrorKind>,
}

#[derive(Clone, Eq, PartialEq, Debug, Fail)]
pub enum VCFParseErrorKind {
    #[fail(display = "I/O Error")]
    IoError,
    #[fail(display = "Header parse error at line: {}  column: {}", line, column)]
    HeaderParseError { line: u32, column: usize },
    #[fail(display = "Record parse error at line: {}  column: {}", line, column)]
    RecordParseError { line: u32, column: usize },
    #[fail(display = "Allele frequency should be number at line: {}", line)]
    FrequencyIsNotNumber { line: u32 },
    #[fail(display = "UTF-8 Error")]
    Utf8Error,
    #[fail(display = "multiple GT record")]
    InvalidGTRecord,
    #[fail(display = "parse number error")]
    ParseNumberError,
    #[fail(display = "parse error")]
    ParseError,
}

impl Fail for VCFParseError {
    fn cause(&self) -> Option<&dyn Fail> {
        self.inner.cause()
    }

    fn backtrace(&self) -> Option<&Backtrace> {
        self.inner.backtrace()
    }
}

impl Display for VCFParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.inner, f)
    }
}

impl VCFParseError {
    pub fn kind(&self) -> VCFParseErrorKind {
        self.inner.get_context().clone()
    }
}

impl From<VCFParseErrorKind> for VCFParseError {
    fn from(kind: VCFParseErrorKind) -> VCFParseError {
        VCFParseError {
            inner: Context::new(kind),
        }
    }
}

impl From<Context<VCFParseErrorKind>> for VCFParseError {
    fn from(inner: Context<VCFParseErrorKind>) -> VCFParseError {
        VCFParseError { inner }
    }
}

impl From<std::io::Error> for VCFParseError {
    fn from(e: std::io::Error) -> VCFParseError {
        e.context(VCFParseErrorKind::IoError).into()
    }
}

impl From<std::str::Utf8Error> for VCFParseError {
    fn from(e: std::str::Utf8Error) -> VCFParseError {
        e.context(VCFParseErrorKind::Utf8Error).into()
    }
}

impl From<std::num::ParseIntError> for VCFParseError {
    fn from(e: std::num::ParseIntError) -> VCFParseError {
        e.context(VCFParseErrorKind::ParseNumberError).into()
    }
}

impl From<std::num::ParseFloatError> for VCFParseError {
    fn from(e: std::num::ParseFloatError) -> VCFParseError {
        e.context(VCFParseErrorKind::ParseNumberError).into()
    }
}

impl From<nom::Err<(&[u8], nom::error::ErrorKind)>> for VCFParseError {
    fn from(_e: nom::Err<(&[u8], nom::error::ErrorKind)>) -> VCFParseError {
        VCFParseErrorKind::ParseError.into()
    }
}

pub fn as_header_error(line: u32, data: &[u8], e: (&[u8], nom::error::ErrorKind)) -> VCFParseError {
    let column = data.offset(e.0);
    VCFParseErrorKind::HeaderParseError { line, column }.into()
}

pub fn as_record_error(line: u32, data: &[u8], e: (&[u8], nom::error::ErrorKind)) -> VCFParseError {
    let column = data.offset(e.0);
    VCFParseErrorKind::RecordParseError { line, column }.into()
}
