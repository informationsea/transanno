use super::Command;
use crate::chain;
use autocompress::{create, open};
use bio::io::fasta::IndexedReader;
use clap::{App, Arg, ArgMatches};
use log::info;

pub struct ChainLeftAlign;

impl Command for ChainLeftAlign {
    fn command_name(&self) -> &'static str {
        "minimap2chain"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Left align and normalize chain file")
            .arg(
                Arg::with_name("original-chain")
                    .index(1)
                    .takes_value(true)
                    .required(true)
                    .help("Original chain file"),
            )
            .arg(
                Arg::with_name("output")
                    .long("output")
                    .short("o")
                    .takes_value(true)
                    .required(true)
                    .help("Output chain file"),
            )
            .arg(
                Arg::with_name("reference_sequence")
                    .long("reference")
                    .short("r")
                    .takes_value(true)
                    .required(true)
                    .help("Reference FASTA (.fai file is required)"),
            )
            .arg(
                Arg::with_name("query_sequence")
                    .long("query")
                    .short("q")
                    .takes_value(true)
                    .required(true)
                    .help("Query FASTA (.fai file is required)"),
            )
    }
    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::LiftOverError> {
        info!("start loading chain");
        let chain_file = open(matches.value_of("original-chain").unwrap())
            .expect("Cannot open input chain file");
        let mut output_file =
            create(matches.value_of("output").unwrap()).expect("Cannot open output chain file");
        let mut reference_seq =
            IndexedReader::from_file(&matches.value_of("reference_sequence").unwrap())
                .expect("Cannot load reference sequence");
        let mut query_seq = IndexedReader::from_file(&matches.value_of("query_sequence").unwrap())
            .expect("Cannot load query sequence");
        let chain_data = chain::ChainFile::load(chain_file).expect("Failed to parse chain file");
        let left_aligned = chain_data
            .left_align(&mut reference_seq, &mut query_seq)
            .expect("Failed to left align");
        left_aligned.write(&mut output_file)?;
        Ok(())
    }
}
