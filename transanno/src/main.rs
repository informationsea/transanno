pub mod commands;

use clap::Parser;
pub use liftover::LiftOverError;
use std::env;

#[derive(Debug, Clone, clap::Parser)]
#[command(version, about = "Transfer annotation to other genome assemblies")]
struct Cli {
    #[command(subcommand)]
    command: commands::Commands,
    #[arg(short = 'v', long = "verbose", action= clap::ArgAction::Count, help="verbose level")]
    verbose: u8,
}

fn main() -> anyhow::Result<()> {
    let app = Cli::parse();

    match app.verbose {
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

    app.command.run()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_clap1() {
        let _cli = Cli::parse_from(&["transanno", "--help"]);
    }

    #[test]
    fn test_clap_liftvcf() {
        let _cli = Cli::parse_from(&["transanno", "liftvcf", "--help"]);
    }

    #[test]
    fn test_clap_liftgene() {
        let _cli = Cli::parse_from(&["transanno", "liftgene", "--help"]);
    }

    #[test]
    fn test_clap_left_align() {
        let _cli = Cli::parse_from(&["transanno", "left-align", "--help"]);
    }

    #[test]
    fn test_clap_minimap2chain() {
        let _cli = Cli::parse_from(&["transanno", "minimap2chain", "--help"]);
    }

    #[test]
    fn test_clap_chain_to_bed_vcf() {
        let _cli = Cli::parse_from(&["transanno", "chain-to-bed-vcf", "--help"]);
    }

    #[test]
    fn test_clap_liftbed() {
        let _cli = Cli::parse_from(&["transanno", "liftbed", "--help"]);
    }
}
