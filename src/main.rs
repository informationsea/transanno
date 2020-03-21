#[macro_use]
extern crate lazy_static;
extern crate failure;

pub mod chain;
mod cli;
pub mod commands;
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
use std::env;

fn main() {
    let app = cli::build_cli();
    let matches = app.get_matches();

    match matches.occurrences_of("verbose") {
        1 => env::set_var("RUST_LOG", "info"),
        2 => env::set_var("RUST_LOG", "debug"),
        3 => env::set_var("RUST_LOG", "trace"),
        _ => {
            if env::var("RUST_LOG").is_err() {
                env::set_var("RUST_LOG", "warn")
            }
        }
    }

    pretty_env_logger::init();

    for one_command in commands::COMMANDS {
        if let Some(matches) = matches.subcommand_matches(one_command.command_name()) {
            one_command.run(matches).expect("Operation Error");
            return;
        }
    }
    unreachable!()
}
