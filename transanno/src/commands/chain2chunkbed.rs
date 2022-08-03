use super::Command;
use crate::cli::validate_integer;
use autocompress::{create, open, CompressionLevel};
use bio::io::fasta::IndexedReader;
use clap::{App, Arg, ArgMatches};
use liftover::chain::{Chain, Chromosome, Strand};
use liftover::{reverse_complement, GenomeSequence};
use log::{debug, trace};
use std::io::{self, Write};
use std::str;

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
                Arg::with_name("reference-vcf")
                    .long("output-reference-vcf")
                    .short("v")
                    .takes_value(true)
                    .help("Output reference/new sequence VCF file (Not sorted)")
            )
            .arg(
                Arg::with_name("query-vcf")
                    .long("output-query-vcf")
                    .short("c")
                    .takes_value(true)
                    .help("Output query/original sequence VCF file (Not sorted)")
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
            ).arg(
                Arg::with_name("svlen")
                    .long("svlen")
                    .short("s")
                    .takes_value(true)
                    .default_value("50")
                    .validator(validate_integer)
                    .help("Do not write nucleotides if a length of reference or alternative sequence is longer than svlen [default: 50]"),
            )
    }
    fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()> {
        chain_to_chunk_bed_helper(
            matches.value_of("chain").unwrap(),
            matches.value_of("reference-sequence").unwrap(),
            matches.value_of("query-sequence").unwrap(),
            matches.value_of("reference-bed").unwrap(),
            matches.value_of("query-bed").unwrap(),
            matches.value_of("reference-vcf"),
            matches.value_of("query-vcf"),
            matches.value_of("svlen").unwrap().parse().unwrap(),
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
    reference_vcf_path: Option<&str>,
    query_vcf_path: Option<&str>,
    sv_len: u64,
) -> anyhow::Result<()> {
    debug!("Reference Sequence: {}", reference_sequence_path);
    debug!("    Query Sequence: {}", query_sequence_path);
    debug!("     Reference BED: {}", reference_bed_path);
    debug!("         Query BED: {}", query_bed_path);
    debug!("     Reference VCF: {:?}", reference_vcf_path);
    debug!("         Query VCF: {:?}", query_vcf_path);

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
    let mut query_vcf_writer = query_vcf_path
        .map(|x| create(x, CompressionLevel::Default).expect("Cannot create Query VCF file"));
    let mut reference_vcf_writer = reference_vcf_path
        .map(|x| create(x, CompressionLevel::Default).expect("Cannot create Reference VCF file"));

    if let Some(query_vcf_writer) = query_vcf_writer.as_mut() {
        write_vcf_header_for_sequence(query_vcf_writer, &query_sequence, query_sequence_path)?;
    }

    if let Some(reference_vcf_writer) = reference_vcf_writer.as_mut() {
        write_vcf_header_for_sequence(
            reference_vcf_writer,
            &reference_sequence,
            reference_sequence_path,
        )?;
    }

    for one_chain in chain_file.chain_list.iter() {
        process_one_chain(
            &one_chain,
            &mut reference_sequence,
            &mut query_sequence,
            &mut reference_bed_writer,
            &mut query_bed_writer,
            &mut reference_vcf_writer,
            &mut query_vcf_writer,
            sv_len,
        )?;
    }

    Ok(())
}

fn process_one_chain<
    W1: Write,
    W2: Write,
    W3: Write,
    W4: Write,
    G1: GenomeSequence,
    G2: GenomeSequence,
>(
    chain: &Chain,
    reference_sequence: &mut G1,
    query_sequence: &mut G2,
    reference_bed_writer: &mut W1,
    query_bed_writer: &mut W2,
    reference_vcf_writer: &mut Option<W3>,
    query_vcf_writer: &mut Option<W4>,
    svlen: u64,
) -> anyhow::Result<()> {
    chain.check_sequence_consistency(reference_sequence, query_sequence)?;
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

    trace!(
        "processing chain / ref: {}:{}-{} ({}/{}/{:?}) / query: {}:{}-{} ({}/{}/{:?})",
        chain.reference_chromosome.name,
        ref_start,
        ref_end,
        chain.reference_strand,
        chain.reference_chromosome.length,
        reference_sequence
            .get_contig_list()
            .iter()
            .filter(|x| x.0 == chain.reference_chromosome.name)
            .next(),
        chain.query_chromosome.name,
        query_start,
        query_end,
        chain.query_strand,
        chain.query_chromosome.length,
        query_sequence
            .get_contig_list()
            .iter()
            .filter(|x| x.0 == chain.query_chromosome.name)
            .next()
    );

    let mut current_reference_pos = chain.reference_start;
    let mut current_query_pos = chain.query_start;
    let mut reference_chunks = Vec::new();
    let mut query_chunks = Vec::new();

    for (i, one_interval) in chain.chain_interval.iter().enumerate() {
        let next_reference_pos = current_reference_pos + one_interval.size;
        let next_query_pos = current_query_pos + one_interval.size;

        let one_ref_chunk = convert_position(
            chain.reference_strand,
            current_reference_pos,
            next_reference_pos,
            chain.reference_chromosome.length,
        );
        reference_chunks.push(one_ref_chunk.clone());

        let one_query_chunk = convert_position(
            chain.query_strand,
            current_query_pos,
            next_query_pos,
            chain.query_chromosome.length,
        );
        query_chunks.push(one_query_chunk.clone());

        trace!(
            "processing chain interval {} / ref: {}:{}-{} / query: {}:{}-{}",
            i,
            chain.reference_chromosome.name,
            one_ref_chunk.0,
            one_ref_chunk.1,
            chain.query_chromosome.name,
            one_query_chunk.0,
            one_query_chunk.1
        );
        process_one_chunk(
            reference_sequence,
            query_sequence,
            reference_vcf_writer,
            query_vcf_writer,
            &chain.reference_chromosome,
            one_ref_chunk.0,
            chain.reference_strand,
            &chain.query_chromosome,
            one_query_chunk.0,
            chain.query_strand,
            one_interval.size,
            svlen,
        )?;

        current_reference_pos = next_reference_pos + one_interval.difference_reference.unwrap_or(0);
        current_query_pos = next_query_pos + one_interval.difference_query.unwrap_or(0);
    }

    reference_chunks.sort();
    query_chunks.sort();

    let ref_thick_start = reference_chunks.iter().fold(String::new(), |mut acc, x| {
        if acc.is_empty() {
            format!("{}", x.0 - ref_start)
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
            format!("{}", x.0 - query_start)
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

    let strand = if chain.reference_strand == chain.query_strand {
        "+"
    } else {
        "-"
    };

    write!(
        reference_bed_writer,
        "{0}\t{1}\t{2}\tchain-{10}-{3}:{4}-{5}\t0\t{6}\t{1}\t{2}\t226,4,27\t{7}\t{8}\t{9}\n",
        chain.reference_chromosome.name,
        ref_start,
        ref_end,
        chain.query_chromosome.name,
        query_start,
        query_end,
        strand,
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
        strand,
        query_chunks.len(),
        query_thick_size,
        query_thick_start,
        chain.chain_id
    )?;

    for (i, (one_query_chunk, one_reference_chunk)) in
        query_chunks.iter().zip(reference_chunks.iter()).enumerate()
    {
        write!(
            reference_bed_writer,
            "{0}\t{1}\t{2}\tchunk-{10}_{11}-{3}:{4}-{5}\t0\t{6}\t{1}\t{2}\t0,123,187\t{7}\t{8}\t{9}\n",
            chain.reference_chromosome.name,
            one_reference_chunk.0,
            one_reference_chunk.1,
            chain.query_chromosome.name,
            one_query_chunk.0,
            one_query_chunk.1,
            strand,
            1,
            one_reference_chunk.1 - one_reference_chunk.0,
            0,
            chain.chain_id,
            i
        )?;

        write!(
            query_bed_writer,
            "{3}\t{4}\t{5}\tchunk-{10}_{11}-{0}:{1}-{2}\t0\t{6}\t{4}\t{5}\t0,123,187\t{7}\t{8}\t{9}\n",
            chain.reference_chromosome.name,
            one_reference_chunk.0,
            one_reference_chunk.1,
            chain.query_chromosome.name,
            one_query_chunk.0,
            one_query_chunk.1,
            strand,
            1,
            one_query_chunk.1 - one_query_chunk.0,
            0,
            chain.chain_id,
            i
        )?;
    }

    Ok(())
}

fn process_one_chunk<G1: GenomeSequence, G2: GenomeSequence, W3: Write, W4: Write>(
    reference_sequence: &mut G1,
    query_sequence: &mut G2,
    reference_vcf_writer: &mut Option<W3>,
    query_vcf_writer: &mut Option<W4>,
    reference_chromosome: &Chromosome,
    reference_start: u64,
    reference_strand: Strand,
    query_chromosome: &Chromosome,
    query_start: u64,
    query_strand: Strand,
    length: u64,
    svlen: u64,
) -> anyhow::Result<()> {
    debug!("Reference Sequence Pointer: {:?}", reference_sequence);
    let reference_seq = reference_sequence.get_sequence(
        &reference_chromosome.name,
        reference_start,
        reference_start + length,
    )?;
    let query_seq =
        query_sequence.get_sequence(&query_chromosome.name, query_start, query_start + length)?;
    let query_seq = if reference_strand == query_strand {
        query_seq
    } else {
        reverse_complement(&query_seq)
    };

    let strand = if reference_strand == query_strand {
        "+"
    } else {
        "-"
    };

    for (i, (ref_base, query_base)) in reference_seq.iter().zip(query_seq.iter()).enumerate() {
        if ref_base != query_base {
            if let Some(reference_vcf_writer) = reference_vcf_writer.as_mut() {
                write!(
                    reference_vcf_writer,
                    "{}\t{}\t.\t{}\t{}\t.\t.\tSTRAND={};TARGET_CRHOM={};TARGET_POS={}\n",
                    reference_chromosome.name,
                    reference_start + (i as u64) + 1,
                    str::from_utf8(&[*ref_base]).unwrap(),
                    str::from_utf8(&[*query_base]).unwrap(),
                    strand,
                    query_chromosome.name,
                    query_start + (i as u64)
                )?;
            }
            if let Some(query_vcf_writer) = query_vcf_writer.as_mut() {
                write!(
                    query_vcf_writer,
                    "{}\t{}\t.\t{}\t{}\t.\t.\tSTRAND={};TARGET_CRHOM={};TARGET_POS={}\n",
                    query_chromosome.name,
                    query_start + (i as u64) + 1,
                    str::from_utf8(&[*query_base]).unwrap(),
                    str::from_utf8(&[*ref_base]).unwrap(),
                    strand,
                    query_chromosome.name,
                    query_start + (i as u64)
                )?;
            }
        }
    }

    Ok(())
}

fn convert_position(strand: Strand, start: u64, end: u64, length: u64) -> (u64, u64) {
    match strand {
        Strand::Forward => (start, end),
        Strand::Reverse => (length - end, length - start),
    }
}

fn write_vcf_header_for_sequence<W: Write, G: GenomeSequence>(
    mut writer: W,
    sequence: &G,
    path: &str,
) -> io::Result<()> {
    let vcf_header = include_bytes!("vcfheader.vcf");
    writer.write_all(vcf_header)?;
    writeln!(writer, "##reference={}", path)?;
    for (one_chrom, chrom_len) in sequence.get_contig_list() {
        writeln!(writer, "##contig=<ID={},length={}>", one_chrom, chrom_len)?;
    }
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

    Ok(())
}
