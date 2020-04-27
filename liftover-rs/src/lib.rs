#[macro_use]
extern crate lazy_static;
extern crate failure;

pub mod chain;
mod defs;
mod error;
pub mod genelift;
pub mod geneparse;
pub mod normalize;
pub mod poslift;
pub mod variantlift;
pub mod vcflift;
pub mod vcfparse;

pub use defs::*;
pub use error::{LiftOverError, LiftOverErrorKind};
