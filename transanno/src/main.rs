mod cli;
pub mod commands;

pub use liftover::LiftOverError;
use log::error;
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
            if let Err(e) = one_command.run(matches) {
                error!("Error: {}", e);
                std::process::exit(1);
            }
            return;
        }
    }
    unreachable!()
}
