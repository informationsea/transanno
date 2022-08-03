use super::Command;
use crate::cli::validate_integer;
use anyhow::Context;
use autocompress::{create, open, CompressionLevel};
use bio::io::fasta::IndexedReader;
use clap::{App, Arg, ArgMatches};
use liftover::{chain, variantlift, vcflift, LiftOverError};
use log::info;

pub struct LiftVcf;

impl Command for LiftVcf {
    fn command_name(&self) -> &'static str {
        "liftvcf"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("LiftOver VCF file")
        .arg(
            Arg::with_name("reference_sequence")
                .long("original-assembly")
                .alias("reference")
                .short("r")
                .takes_value(true)
                .required(true)
                .help("Original assembly FASTA (.fai file is required)"),
        )
        .arg(
            Arg::with_name("query_sequence")
                .long("new-assembly")
                .alias("query")
                .short("q")
                .takes_value(true)
                .required(true)
                .help("New assembly FASTA (.fai file is required)"),
        )
        .arg(
            Arg::with_name("chain")
                .long("chain")
                .short("c")
                .takes_value(true)
                .required(true)
                .help("chain file"),
        )
        .arg(
            Arg::with_name("vcf")
                .long("vcf")
                .short("f")
                .takes_value(true)
                .required(true)
                .help("input VCF file to liftOver"),
        )
        .arg(
            Arg::with_name("output")
                .long("output")
                .short("o")
                .takes_value(true)
                .required(true)
                .help("output VCF file for succeeded to liftOver records (This file is not sorted)"),
        )
        .arg(
            Arg::with_name("fail")
                .long("fail")
                .short("a")
                .takes_value(true)
                .required(true)
                .help("output VCF file for failed to liftOver records"),
        ).arg(
            Arg::with_name("allow-multimap")
                .long("allow-multi-map")
                .short("m")
                .help("Allow multi-map")
        ).arg(
            Arg::with_name("acceptable-deletion")
                .long("acceptable-deletion")
                .short("d")
                .default_value("3")
                .takes_value(true)
                .validator(validate_integer)
                .help("length of acceptable deletion")
        ).arg(
            Arg::with_name("acceptable-insertion")
                .long("acceptable-insertion")
                .short("i")
                .default_value("3")
                .takes_value(true)
                .validator(validate_integer)
                .help("length of acceptable insertion")
        ).arg(
            Arg::with_name("do-not-rewrite-info")
                .long("no-rewrite-info")
                .help("Do not rewrite order of INFO tags")
        ).arg(
            Arg::with_name("do-not-rewrite-format")
                .long("no-rewrite-format")
                .help("Do not rewrite order of FORMAT tags")
        ).arg(
            Arg::with_name("do-not-rewrite-gt")
                .long("no-rewrite-gt")
                .help("Do not rewrite order of GT")
        ).arg(
            Arg::with_name("do-not-rewrite-allele-frequency")
                .long("no-rewrite-allele-frequency")
                .help("Do not rewrite AF or other allele frequency info")
        ).arg(
            Arg::with_name("do-not-rewrite-allele-count")
                .long("no-rewrite-allele-count")
                .help("Do not rewrite AC or other count frequency info")
        ).arg(
            Arg::with_name("do-not-swap-ref-alt")
                .long("noswap")
                .help("Do not swap ref/alt when reference allele is changed. This option is suitable to do liftOver clinVar, COSMIC annotations")
        ).arg(
            Arg::with_name("do-not-left-align-chain")
                .long("no-left-align-chain")
                .help("Do not run left align chain file")
        ).arg(
            Arg::with_name("do_not_use_dot_when_alt_equal_to_ref")
                .long("do-not-use-dot-when-alt-equal-to-ref")
                .help("Do not use dot as ALT when ALT column is equal to REF")
        ).arg(
            Arg::with_name("do_not_prefer_cis_contig_when_multimap")
                .long("do-not-prefer-same-contig-when-multimap")
                .help("Do not prefer same name contig when a variant lifted into multiple positions. (When you use this option, a variant which lifted into a main chromosome and alternative contigs, lift over will be failed if multimap is not allowed)")
        )
    }

    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        info!("start loading chain and fasta");
        let mut reference_seq =
            IndexedReader::from_file(&matches.value_of("reference_sequence").unwrap())
                .context("Failed to load reference sequence")?;
        let mut query_seq = IndexedReader::from_file(&matches.value_of("query_sequence").unwrap())
            .context("Failed to load query sequence")?;
        let chain =
            chain::ChainFile::load(autocompress::open(matches.value_of("chain").unwrap())?)?
                .left_align(&mut reference_seq, &mut query_seq)
                .context("Failed to load chain file")?;
        // Reference/Query sequence and chain consistency
        for one_chain in chain.chain_list.iter() {
            match one_chain.check_sequence_consistency(&mut reference_seq, &mut query_seq) {
                Ok(_) => (),
                Err(LiftOverError::ChromosomeNotFound(_)) => (),
                Err(e) => return Err(e.into()),
            }
        }

        let variant_liftover = variantlift::VariantLiftOver::new(chain, reference_seq, query_seq);
        let mut vcf_lift = vcflift::VCFLiftOver::new(
            variant_liftover,
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
                .do_not_rewrite_allele_frequency(
                    matches.is_present("do-not-rewrite-allele-frequency"),
                )
                .do_not_rewrite_allele_count(matches.is_present("do-not-rewrite-allele-count"))
                .do_not_swap_ref_alt(matches.is_present("do-not-swap-ref-alt"))
                .do_not_left_align_chain_file(matches.is_present("do-not-left-align-chain"))
                .do_not_use_dot_when_alt_equal_to_ref(
                    matches.is_present("do_not_use_dot_when_alt_equal_to_ref"),
                )
                .do_not_prefer_cis_contig_when_multimap(
                    matches.is_present("do_not_prefer_cis_contig_when_multimap"),
                ),
        );
        info!("chain file and fasta files were loaded");

        let uncompressed_reader =
            open(matches.value_of("vcf").unwrap()).expect("Cannot open input VCF");
        let success_writer = create(
            matches.value_of("output").unwrap(),
            CompressionLevel::Default,
        )
        .expect("Cannot create output VCF");
        let failed_writer = create(matches.value_of("fail").unwrap(), CompressionLevel::Default)
            .expect("Cannot create output VCF for failed records");
        vcf_lift.lift_vcf(uncompressed_reader, success_writer, failed_writer)?;
        Ok(())
    }
}
