use failure::*;
use std::fmt::{self, Display};

#[derive(Debug)]
pub struct LiftOverError {
    inner: failure::Context<LiftOverErrorKind>,
}

#[derive(Clone, Eq, PartialEq, Debug, Fail)]
pub enum LiftOverErrorKind {
    #[fail(display = "I/O Error")]
    IoError,
    #[fail(display = "UTF-8 or IO Error")]
    Utf8IoError,
    #[fail(display = "VCF parse error")]
    VCFParseError,
    #[fail(display = "Invalid number of header at line {}", _0)]
    InvalidNumberOfHeader(u32),
    #[fail(
        display = "Invalid header at line {}. A header line should starts with \"chain\"",
        _0
    )]
    NoChainHeaderFound(u32),
    #[fail(display = "Invalid strand at line {}", _0)]
    InvalidStrand(u32),
    #[fail(display = "Invalid chromosome length at line {}", _0)]
    InvalidChromosomeLength(u32),
    #[fail(display = "Invalid number of columns at line {}", _0)]
    InvalidNumberOfColumns(u32),
    #[fail(display = "Reference sequence is different from expected sequence")]
    DifferentReference,
    #[fail(display = "Parse integer error")]
    ParseIntError,
    #[fail(display = "Parse strand error")]
    ParseStrandError,
    #[fail(display = "stand for reference should be forward")]
    ReferenceStrandShouldForward,
    #[fail(
        display = "length of chromosome {} is not equal to length in chain file. Are you using correct reference?",
        _0
    )]
    ReferenceChromosomeLengthIsNotMatch(String),
    #[fail(
        display = "length of chromosome {} is not equal to length in chain file. Are you using correct query?",
        _0
    )]
    QueryChromosomeLengthIsNotMatch(String),
    #[fail(display = "unknown sequence error: {}", _0)]
    UnknownSequenceError(String),
    #[fail(display = "Failed to parse gene annotation")]
    GeneParseError,
}

impl Fail for LiftOverError {
    fn cause(&self) -> Option<&dyn Fail> {
        self.inner.cause()
    }

    fn backtrace(&self) -> Option<&Backtrace> {
        self.inner.backtrace()
    }
}

impl Display for LiftOverError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Display::fmt(&self.inner, f)
    }
}

impl LiftOverError {
    pub fn kind(&self) -> LiftOverErrorKind {
        self.inner.get_context().clone()
    }
}

impl From<LiftOverErrorKind> for LiftOverError {
    fn from(kind: LiftOverErrorKind) -> LiftOverError {
        LiftOverError {
            inner: Context::new(kind),
        }
    }
}

impl From<Context<LiftOverErrorKind>> for LiftOverError {
    fn from(inner: Context<LiftOverErrorKind>) -> LiftOverError {
        LiftOverError { inner }
    }
}

impl From<std::io::Error> for LiftOverError {
    fn from(e: std::io::Error) -> LiftOverError {
        e.context(LiftOverErrorKind::IoError).into()
    }
}

impl From<csv::Error> for LiftOverError {
    fn from(e: csv::Error) -> LiftOverError {
        e.context(LiftOverErrorKind::Utf8IoError).into()
    }
}

impl From<std::num::ParseIntError> for LiftOverError {
    fn from(e: std::num::ParseIntError) -> LiftOverError {
        e.context(LiftOverErrorKind::ParseIntError).into()
    }
}

impl From<crate::vcfparse::VCFParseError> for LiftOverError {
    fn from(e: crate::vcfparse::VCFParseError) -> LiftOverError {
        e.context(LiftOverErrorKind::IoError).into()
    }
}

impl From<crate::geneparse::GeneParseError> for LiftOverError {
    fn from(e: crate::geneparse::GeneParseError) -> LiftOverError {
        e.context(LiftOverErrorKind::GeneParseError).into()
    }
}
