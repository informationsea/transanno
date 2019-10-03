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

use bio::io::fasta::IndexedReader;
use clap::{crate_name, ArgMatches, Shell};
pub use defs::*;
pub use error::{LiftOverError, LiftOverErrorKind};
use log::info;
use std::env;
use std::fs;
use std::io::BufRead;

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

    if let Some(matches) = matches.subcommand_matches("generate-completions") {
        generate_completions(matches);
    } else if let Some(matches) = matches.subcommand_matches("liftvcf") {
        lift_vcf(matches);
    } else if let Some(matches) = matches.subcommand_matches("chain-left-align") {
        left_align_chain(matches);
    } else if let Some(matches) = matches.subcommand_matches("chain-to-bed-vcf") {
        commands::chain_to_bed_vcf(matches)
    } else if let Some(matches) = matches.subcommand_matches("minimap2-to-chain") {
        commands::minimap2_to_chain(matches)
    } else if let Some(matches) = matches.subcommand_matches("lift-gff") {
        commands::lift_gene(matches)
    } else {
        unreachable!()
    }
}

fn left_align_chain(matches: &ArgMatches) {
    info!("start loading chain");
    let chain_file = adaptive_open(matches.value_of("original-chain").unwrap())
        .expect("Cannot open input chain file");
    let mut output_file = adaptive_create(matches.value_of("output").unwrap())
        .expect("Cannot open output chain file");
    let mut reference_seq =
        IndexedReader::from_file(&matches.value_of("reference_sequence").unwrap())
            .expect("Cannot load reference sequence");
    let mut query_seq = IndexedReader::from_file(&matches.value_of("query_sequence").unwrap())
        .expect("Cannot load query sequence");
    let chain_data = chain::ChainFile::load(chain_file).expect("Failed to parse chain file");
    let left_aligned = chain_data
        .left_align(&mut reference_seq, &mut query_seq)
        .expect("Failed to left align");
    left_aligned
        .write(&mut output_file)
        .expect("Failed to write left aligned chain");
}

fn lift_vcf(matches: &ArgMatches) {
    info!("start loading chain and fasta");
    let mut vcf_lift = vcflift::VCFLiftOver::load(
        matches.value_of("chain").unwrap(),
        matches.value_of("reference_sequence").unwrap(),
        matches.value_of("query_sequence").unwrap(),
        vcflift::VCFLiftOverParameters::new()
            .allow_multimap(matches.is_present("allow-multimap"))
            .acceptable_deletion(
                matches
                    .value_of("acceptable-deletion")
                    .unwrap()
                    .parse()
                    .expect("length of acceptable deletion should be larger than zero"),
            )
            .acceptable_insertion(
                matches
                    .value_of("acceptable-insertion")
                    .unwrap()
                    .parse()
                    .expect("length of acceptable deletion should be larger than zero"),
            )
            .do_not_rewrite_info(matches.is_present("do-not-rewrite-info"))
            .do_not_rewrite_format(matches.is_present("do-not-rewrite-format"))
            .do_not_rewrite_gt(matches.is_present("do-not-rewrite-gt"))
            .do_not_rewrite_allele_frequency(matches.is_present("do-not-rewrite-allele-frequency"))
            .do_not_rewrite_allele_count(matches.is_present("do-not-rewrite-allele-count"))
            .do_not_swap_ref_alt(matches.is_present("do-not-swap-ref-alt"))
            .do_not_left_align_chain_file(matches.is_present("do-not-left-align-chain"))
            .do_not_use_dot_when_alt_equal_to_ref(
                matches.is_present("do_not_use_dot_when_alt_equal_to_ref"),
            )
            .do_not_prefer_cis_contig_when_multimap(
                matches.is_present("do_not_prefer_cis_contig_when_multimap"),
            ),
    )
    .expect("Cannot load chain/FASTA file");
    info!("chain file and fasta files were loaded");

    let uncompressed_reader: Box<dyn BufRead> =
        adaptive_open(matches.value_of("vcf").unwrap()).expect("Cannot open input VCF");
    let success_writer =
        adaptive_create(matches.value_of("output").unwrap()).expect("Cannot create output VCF");
    let failed_writer = adaptive_create(matches.value_of("fail").unwrap())
        .expect("Cannot create output VCF for failed records");
    vcf_lift
        .lift_vcf(uncompressed_reader, success_writer, failed_writer)
        .expect("Failed to lift over");
}

fn generate_completions(matches: &ArgMatches) {
    let outdir = matches.value_of("output").unwrap();
    fs::create_dir_all(&outdir).expect("Cannot create output directory");
    let mut app = cli::build_cli();
    app.gen_completions(crate_name!(), Shell::Bash, &outdir);
    app.gen_completions(crate_name!(), Shell::Zsh, &outdir);
    app.gen_completions(crate_name!(), Shell::PowerShell, &outdir);
    app.gen_completions(crate_name!(), Shell::Fish, &outdir);
}
