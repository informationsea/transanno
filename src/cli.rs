use clap::{crate_authors, crate_version, App, Arg, SubCommand};

include!(concat!(env!("OUT_DIR"), "/git_version.rs"));

pub fn build_cli() -> App<'static, 'static> {
    App::new("transanno")
        .version(crate_version!())
        .long_version(concat!(crate_version!(), " git:", git_version!()))
        .author(crate_authors!())
        .about("Transfer annotation to other genome assemblies")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Verbose mode (-v: info, -vv: debug, -vvv: trace)"),
        )
        .subcommand(
            SubCommand::with_name("minimap2-to-chain")
            .version(crate_version!())
            .about("Convert minimap2 result to chain file")
            .long_about(r#"Convert minimap2 result to chain file

A paf file should be created with a command shown in below.

$ minimap2 -cx asm5 --cs QUERY_FASTA REFERENCE_FASTA > PAF_FILE.paf
"#)
            .arg(
                Arg::with_name("paf")
                    .index(1)
                    .takes_value(true)
                    .required(true)
                    .help("Input paf file")
            )
            .arg(
                Arg::with_name("chain")
                    .long("output")
                    .short("o")
                    .takes_value(true)
                    .required(true)
                    .help("Output chain file")
            )
        )
        .subcommand(
            SubCommand::with_name("chain-to-bed-vcf")
            .version(crate_version!())
            .about("Create BED and VCF file from chain file")
            .arg(
                Arg::with_name("chain")
                    .index(1)
                    .takes_value(true)
                    .required(true)
                    .help("Input Chain file")
            )
            .arg(
                Arg::with_name("reference-bed")
                    .long("output-reference-bed")
                    .short("b")
                    .takes_value(true)
                    .required(true)
                    .help("Output Reference BED file (Not sorted)")
            )
            .arg(
                Arg::with_name("query-bed")
                    .long("output-query-bed")
                    .short("d")
                    .takes_value(true)
                    .required(true)
                    .help("Output Query BED file (Not sorted)")
            )
            .arg(
                Arg::with_name("reference-vcf")
                    .long("output-reference-vcf")
                    .short("v")
                    .takes_value(true)
                    .required(true)
                    .help("Output Reference VCF file (Not sorted)")
            )
            .arg(
                Arg::with_name("query-vcf")
                    .long("output-query-vcf")
                    .short("c")
                    .takes_value(true)
                    .required(true)
                    .help("Output Query VCF file (Not sorted)")
            )
            .arg(
                Arg::with_name("reference-sequence")
                    .long("reference")
                    .short("r")
                    .takes_value(true)
                    .required(true)
                    .help("Reference FASTA (.fai file is required)"),
            )
            .arg(
                Arg::with_name("query-sequence")
                    .long("query")
                    .short("q")
                    .takes_value(true)
                    .required(true)
                    .help("Query FASTA (.fai file is required)"),
            ).arg(
                Arg::with_name("svlen")
                    .long("svlen")
                    .short("s")
                    .takes_value(true)
                    .default_value("50")
                    .validator(validate_integer)
                    .help("Do not write nucleotides if a length of reference or alternative sequence is longer than svlen [default: 50]"),
            )
        )
        .subcommand(
            SubCommand::with_name("chain-left-align")
                .version(crate_version!())
                .about("Left align and normalize chain file")
                .arg(
                    Arg::with_name("original-chain")
                        .index(1)
                        .takes_value(true)
                        .required(true)
                        .help("Original chain file")
                )
                .arg(
                    Arg::with_name("output")
                        .long("output")
                        .short("o")
                        .takes_value(true)
                        .required(true)
                        .help("Output chain file")
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

        )
        .subcommand(
            SubCommand::with_name("liftvcf")
                .version(crate_version!())
                .about("LiftOver VCF file")
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
                ),
        )
        .subcommand(
            SubCommand::with_name("liftgene").version(crate_version!())
                .about("Lift GENCODE or Ensemble GFF3/GTF file")
                .arg(
                    Arg::with_name("chain")
                        .long("chain")
                        .short("c")
                        .required(true)
                        .takes_value(true)
                        .help("chain file")
                ).arg(
                    Arg::with_name("gff")
                        .index(1)
                        .required(true)
                        .takes_value(true)
                        .help("input GFF3/GTF file (GENCODE/Ensemble)")
                ).arg(
                    Arg::with_name("format")
                        .long("format")
                        .possible_values(&["auto", "GFF3", "GTF"])
                        .takes_value(true)
                        .default_value("auto")
                        .help("Input file format")
                ).arg(
                    Arg::with_name("output")
                        .long("output")
                        .short("o")
                        .required(true)
                        .takes_value(true)
                        .help("GFF3/GTF output path (unsorted)")
                ).arg(
                    Arg::with_name("failed")
                        .long("failed")
                        .short("f")
                        .required(true)
                        .takes_value(true)
                        .help("Failed to liftOver GFF3/GTF output path")
                )
        )
        .subcommand(
            SubCommand::with_name("generate-completions").version(crate_version!())
                .about("Generate completion files")
                .arg(
                    Arg::with_name("output")
                        .index(1)
                        .help("Output directory")
                        .required(true)
                        .takes_value(true)
                )
        )
        .setting(clap::AppSettings::SubcommandRequired)
        .global_setting(clap::AppSettings::ColoredHelp)
        .global_setting(clap::AppSettings::ColorAuto)
}

fn validate_integer(text: String) -> Result<(), String> {
    for one in text.chars() {
        if !one.is_ascii_digit() {
            return Err(format!("{} is not number", text));
        }
    }
    Ok(())
}
