use failure::*;

#[derive(Clone, Eq, PartialEq, Debug, Fail, PartialOrd, Ord)]
pub enum VariantLiftOverError {
    #[fail(display = "Unacceptable large deletion")]
    UnacceptableLargeDeletion {
        chromosome: String,
        start: u64,
        end: u64,
    },
    #[fail(display = "Unacceptable large insertion")]
    UnacceptableLargeInsertion {
        chromosome: String,
        start: u64,
        end: u64,
    },
    #[fail(display = "Unknown sequence name: {}", _0)]
    UnknownSequenceName(String),
    #[fail(display = "REF is different from expected reference sequence")]
    ReferenceSequenceIsNotMatch,
}
