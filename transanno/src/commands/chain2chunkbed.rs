use super::Command;
use crate::cli::validate_integer;
use autocompress::{create, open, CompressionLevel};
use bio::io::fasta::IndexedReader;
use clap::{App, Arg, ArgMatches};
use liftover::chain::{Chain, Strand};
use liftover::LiftOverError;
use liftover::{reverse_acid, reverse_complement, GenomeSequence};
use std::io::{self, BufRead, Read, Seek, Write};

pub struct Chain2ChunkBed;

impl Command for Chain2ChunkBed {
    fn command_name(&self) -> &'static str {
        "chain-to-chunk-bed"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Create BED and VCF file from chain file")
            .arg(
                Arg::with_name("chain")
                    .index(1)
                    .takes_value(true)
                    .required(true)
                    .help("Input Chain file"),
            )
            .arg(
                Arg::with_name("reference-bed")
                    .long("output-reference-bed")
                    .short("b")
                    .takes_value(true)
                    .required(true)
                    .help("Output reference/new sequence BED file (Not sorted)"),
            )
            .arg(
                Arg::with_name("query-bed")
                    .long("output-query-bed")
                    .short("d")
                    .takes_value(true)
                    .required(true)
                    .help("Output query/original sequence BED file (Not sorted)"),
            )
            .arg(
                Arg::with_name("reference-sequence")
                    .long("reference")
                    .short("r")
                    .takes_value(true)
                    .required(true)
                    .help("Reference/new sequence FASTA (.fai file is required)"),
            )
            .arg(
                Arg::with_name("query-sequence")
                    .long("query")
                    .short("q")
                    .takes_value(true)
                    .required(true)
                    .help("Query/original sequence FASTA (.fai file is required)"),
            )
    }
    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        chain_to_chunk_bed_helper(
            matches.value_of("chain").unwrap(),
            matches.value_of("reference-sequence").unwrap(),
            matches.value_of("query-sequence").unwrap(),
            matches.value_of("reference-bed").unwrap(),
            matches.value_of("query-bed").unwrap(),
        )?;
        Ok(())
    }
}

fn chain_to_chunk_bed_helper(
    chain_path: &str,
    reference_sequence_path: &str,
    query_sequence_path: &str,
    reference_bed_path: &str,
    query_bed_path: &str,
) -> anyhow::Result<()> {
    let mut reference_sequence =
        IndexedReader::from_file(&reference_sequence_path).expect("Cannot open reference sequence");
    let mut query_sequence =
        IndexedReader::from_file(&query_sequence_path).expect("Cannot open query sequence");
    let chain_file =
        liftover::chain::ChainFile::load(open(chain_path).expect("Cannot open chain file"))
            .expect("Cannot parse chain file");
    let mut query_bed_writer =
        create(query_bed_path, CompressionLevel::Default).expect("Cannot create Query BED file");
    let mut reference_bed_writer = create(reference_bed_path, CompressionLevel::Default)
        .expect("Cannot create Reference BED file");

    for (i, one_chain) in chain_file.chain_list.iter().enumerate() {
        one_chain_to_bed(
            &one_chain,
            &mut reference_sequence,
            &mut query_sequence,
            &mut reference_bed_writer,
            &mut query_bed_writer,
        )?;
    }

    Ok(())
}

fn one_chain_to_bed<W1: Write, W2: Write, R1: Read + Seek, R2: Read + Seek>(
    chain: &Chain,
    reference_sequence: &mut IndexedReader<R1>,
    query_sequence: &mut IndexedReader<R2>,
    reference_bed_writer: &mut W1,
    query_bed_writer: &mut W2,
) -> anyhow::Result<()> {
    let (ref_start, ref_end) = convert_position(
        chain.reference_strand,
        chain.reference_start,
        chain.reference_end,
        chain.reference_chromosome.length,
    );
    let (query_start, query_end) = convert_position(
        chain.query_strand,
        chain.query_start,
        chain.query_end,
        chain.query_chromosome.length,
    );

    let mut current_reference_pos = chain.reference_start;
    let mut current_query_pos = chain.query_start;
    let mut reference_chunks = Vec::new();
    let mut query_chunks = Vec::new();
    for one_interval in chain.chain_interval.iter() {
        let next_reference_pos = current_reference_pos + one_interval.size;
        let next_query_pos = current_query_pos + one_interval.size;

        let one_ref_chunk = convert_position(
            chain.reference_strand,
            current_reference_pos,
            next_reference_pos,
            chain.reference_chromosome.length,
        );
        reference_chunks.push(one_ref_chunk);

        let one_query_chunk = convert_position(
            chain.query_strand,
            current_query_pos,
            next_query_pos,
            chain.query_chromosome.length,
        );
        query_chunks.push(one_query_chunk);

        current_reference_pos = next_reference_pos + one_interval.difference_reference.unwrap_or(0);
        current_query_pos = next_query_pos + one_interval.difference_query.unwrap_or(0);
    }

    reference_chunks.sort();
    query_chunks.sort();

    let ref_thick_start = reference_chunks.iter().fold(String::new(), |mut acc, x| {
        if acc.is_empty() {
            format!("{}", x.0)
        } else {
            acc.push_str(&format!(",{}", x.0 - ref_start));
            acc
        }
    });
    let ref_thick_size = reference_chunks.iter().fold(String::new(), |mut acc, x| {
        if acc.is_empty() {
            format!("{}", x.1 - x.0)
        } else {
            acc.push_str(&format!(",{}", x.1 - x.0));
            acc
        }
    });
    let query_thick_start = query_chunks.iter().fold(String::new(), |mut acc, x| {
        if acc.is_empty() {
            format!("{}", x.0)
        } else {
            acc.push_str(&format!(",{}", x.0 - query_start));
            acc
        }
    });
    let query_thick_size = query_chunks.iter().fold(String::new(), |mut acc, x| {
        if acc.is_empty() {
            format!("{}", x.1 - x.0)
        } else {
            acc.push_str(&format!(",{}", x.1 - x.0));
            acc
        }
    });

    write!(
        reference_bed_writer,
        "{0}\t{1}\t{2}\tchain-{10}-{3}:{4}-{5}\t0\t{6}\t{1}\t{2}\t226,4,27\t{7}\t{8}\t{9}\n",
        chain.reference_chromosome.name,
        ref_start,
        ref_end,
        chain.query_chromosome.name,
        query_start,
        query_end,
        if chain.reference_strand == chain.query_strand {
            "+"
        } else {
            "-"
        },
        reference_chunks.len(),
        ref_thick_size,
        ref_thick_start,
        chain.chain_id
    )?;

    write!(
        query_bed_writer,
        "{3}\t{4}\t{5}\tchain-{10}-{0}:{1}-{2}\t0\t{6}\t{4}\t{5}\t226,4,27\t{7}\t{8}\t{9}\n",
        chain.reference_chromosome.name,
        ref_start,
        ref_end,
        chain.query_chromosome.name,
        query_start,
        query_end,
        if chain.reference_strand == chain.query_strand {
            "+"
        } else {
            "-"
        },
        query_chunks.len(),
        query_thick_size,
        query_thick_start,
        chain.chain_id
    )?;

    Ok(())
}

fn convert_position(strand: Strand, start: u64, end: u64, length: u64) -> (u64, u64) {
    match strand {
        Strand::Forward => (start, end),
        Strand::Reverse => (length - end, length - start),
    }
}
