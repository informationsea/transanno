use thiserror::Error;

#[derive(Debug, Error)]
pub enum VariantLiftOverError {
    #[error("Unacceptable large deletion")]
    UnacceptableLargeDeletion {
        chromosome: String,
        start: u64,
        end: u64,
    },
    #[error("Unacceptable large insertion")]
    UnacceptableLargeInsertion {
        chromosome: String,
        start: u64,
        end: u64,
    },
    #[error("Unknown sequence name: {}", _0)]
    UnknownSequenceName(String),
    #[error("REF is different from expected reference sequence")]
    ReferenceSequenceIsNotMatch,
}
