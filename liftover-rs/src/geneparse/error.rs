use failure::*;
use log::debug;
use nom::error::ErrorKind;
use nom::Err;
use std::fmt::{self, Display};
use std::num::ParseIntError;

#[derive(Debug)]
pub struct GeneParseError {
    inner: failure::Context<GeneParseErrorKind>,
}

#[derive(Clone, Eq, PartialEq, Debug, Fail)]
pub enum GeneParseErrorKind {
    #[fail(display = "I/O Error")]
    IoError,
    #[fail(display = "Parse error")]
    ParseError,
    #[fail(display = "Parse error at line: {}", _0)]
    ParseErrorAtLine(u64),
    #[fail(display = "Gene grouping error")]
    GroupingError,
    #[fail(display = "# of refGene columns should be 16")]
    RefGeneInvalidColumnNumber,
    #[fail(display = "Invalid strand")]
    StrandParseError,
    #[fail(display = "CDS stat parse error")]
    CDSStatParseError,
}

impl Fail for GeneParseError {
    fn cause(&self) -> Option<&dyn Fail> {
        self.inner.cause()
    }

    fn backtrace(&self) -> Option<&Backtrace> {
        self.inner.backtrace()
    }
}

impl Display for GeneParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.inner, f)
    }
}

impl GeneParseError {
    pub fn kind(&self) -> GeneParseErrorKind {
        self.inner.get_context().clone()
    }
}

impl From<GeneParseErrorKind> for GeneParseError {
    fn from(kind: GeneParseErrorKind) -> GeneParseError {
        GeneParseError {
            inner: Context::new(kind),
        }
    }
}

impl From<Context<GeneParseErrorKind>> for GeneParseError {
    fn from(inner: Context<GeneParseErrorKind>) -> GeneParseError {
        GeneParseError { inner }
    }
}

impl From<std::io::Error> for GeneParseError {
    fn from(e: std::io::Error) -> GeneParseError {
        e.context(GeneParseErrorKind::IoError).into()
    }
}

impl From<Err<(&str, ErrorKind)>> for GeneParseError {
    fn from(e: Err<(&str, ErrorKind)>) -> GeneParseError {
        println!("nom error: {:?}", e);
        GeneParseErrorKind::ParseError.into()
    }
}

impl From<ParseIntError> for GeneParseError {
    fn from(e: ParseIntError) -> GeneParseError {
        debug!("parse int error: {:?}", e);
        GeneParseErrorKind::ParseError.into()
    }
}
