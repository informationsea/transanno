use crate::utils::{create, open};
use anyhow::Context;
use bio::io::fasta::IndexedReader;
use clap::Args;
use liftover::chain::{Chain, Strand};
use liftover::LiftOverError;
use liftover::{reverse_acid, reverse_complement, GenomeSequence};
use log::warn;
use std::io;
use std::io::Write;

#[derive(Debug, Clone, Args)]
#[command(about = "Create BED and VCF file from chain file")]
pub struct Chain2BedVcf {
    #[arg(help = "Input Chain file")]
    chain: String,
    #[arg(
        long = "output-original-bed",
        alias = "output-reference-bed",
        short = 'b',
        help = "Output original assembly BED file (Not sorted)"
    )]
    original_bed: String,
    #[arg(
        long = "output-new-bed",
        alias = "output-query-bed",
        short = 'd',
        help = "Output new assembly BED file (Not sorted)"
    )]
    new_bed: String,
    #[arg(
        long = "output-original-vcf",
        alias = "output-reference-vcf",
        short = 'v',
        help = "Output original assembly VCF file (Not sorted)"
    )]
    original_vcf: String,
    #[arg(
        long = "output-new-vcf",
        alias = "output-query-vcf",
        short = 'c',
        help = "Output new assembly VCF file (Not sorted)"
    )]
    new_vcf: String,
    #[arg(
        long = "original",
        alias = "reference",
        short = 'r',
        help = "Original assembly FASTA (.fai file is required)"
    )]
    original_sequence: String,
    #[arg(
        long = "new",
        alias = "query",
        short = 'q',
        help = "New assembly FASTA (.fai file is required)"
    )]
    new_sequence: String,
    #[arg(
        long = "svlen",
        short = 's',
        default_value = "50",
        help = "Do not write nucleotides if a length of reference or alternative sequence is longer than svlen"
    )]
    svlen: usize,
}

impl Chain2BedVcf {
    pub fn run(&self) -> anyhow::Result<()> {
        chain_to_bed_vcf_helper(
            &self.chain,
            &self.original_sequence,
            &self.new_sequence,
            &self.original_vcf,
            &self.new_vcf,
            &self.original_bed,
            &self.new_bed,
            self.svlen,
        )?;
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn chain_to_bed_vcf_helper(
    chain_path: &str,
    original_sequence_path: &str,
    new_sequence_path: &str,
    original_vcf_path: &str,
    new_vcf_path: &str,
    original_bed_path: &str,
    new_bed_path: &str,
    sv_len: usize,
) -> anyhow::Result<()> {
    let mut original_sequence = IndexedReader::from_file(&original_sequence_path)
        .with_context(|| format!("Cannot open original sequence: {original_sequence_path}"))?;
    let mut new_sequence = IndexedReader::from_file(&new_sequence_path)
        .with_context(|| format!("Cannot open query sequence: {new_sequence_path}"))?;
    let chain_file = liftover::chain::ChainFile::load(
        open(chain_path).with_context(|| format!("Cannot open chain file: {chain_path}"))?,
    )
    .with_context(|| format!("Cannot parse chain file: {chain_path}"))?;

    let mut new_vcf_writer = create(new_vcf_path)
        .with_context(|| format!("Cannot create new assembly VCF file: {new_vcf_path}"))?;
    write_vcf_header_for_sequence(&mut new_vcf_writer, &new_sequence, new_sequence_path)
        .context("Cannot write new assembly VCF")?;
    let mut new_bed_writer = create(new_bed_path)
        .with_context(|| format!("Cannot create new assembly BED file: {new_bed_path}"))?;

    let mut original_vcf_writer = create(original_vcf_path).with_context(|| {
        format!("Cannot create original assembly VCF file: {original_vcf_path}")
    })?;
    write_vcf_header_for_sequence(
        &mut original_vcf_writer,
        &original_sequence,
        original_sequence_path,
    )
    .context("Cannot write original assembly VCF")?;

    let mut original_bed_writer = create(original_bed_path).with_context(|| {
        format!("Cannot create original assembly BED file: {original_bed_path}")
    })?;

    for one_chain in chain_file.chain_list {
        // Check chromosome length
        match one_chain.check_sequence_consistency(&mut original_sequence, &mut new_sequence) {
            Ok(_) => (),
            Err(LiftOverError::ChromosomeNotFound(_)) => {
                warn!("Skip chain ID: {}", one_chain.chain_id);
                continue;
            }
            Err(e) => return Err(e.into()),
        }

        // write BED
        let (reference_start, reference_end) = convert_position(
            one_chain.original_start,
            one_chain.original_end,
            one_chain.original_chromosome.length,
            one_chain.original_strand,
        );
        let (query_start, query_end) = convert_position(
            one_chain.new_start,
            one_chain.new_end,
            one_chain.new_chromosome.length,
            one_chain.new_strand,
        );

        writeln!(
            original_bed_writer,
            "{}\t{}\t{}\tchain_id:{};{}:{}-{}\t{}\t{}",
            one_chain.original_chromosome.name,
            reference_start,
            reference_end,
            one_chain.chain_id,
            one_chain.new_chromosome.name,
            query_start + 1,
            query_end,
            one_chain.score,
            if one_chain.original_strand == one_chain.new_strand {
                "+"
            } else {
                "-"
            },
        )
        .expect("Cannot write reference BED");

        writeln!(
            new_bed_writer,
            "{}\t{}\t{}\tchain_id:{};{}:{}-{}\t{}\t{}",
            one_chain.new_chromosome.name,
            query_start,
            query_end,
            one_chain.chain_id,
            one_chain.original_chromosome.name,
            reference_start + 1,
            reference_end,
            one_chain.score,
            if one_chain.original_strand == one_chain.new_strand {
                "+"
            } else {
                "-"
            },
        )
        .expect("Cannot write query BED");

        write_vcf_entry(
            one_chain,
            &mut original_sequence,
            &mut new_sequence,
            &mut original_vcf_writer,
            &mut new_vcf_writer,
            sv_len,
        )?;
    }

    Ok(())
}

fn write_vcf_entry<G: GenomeSequence, W: io::Write>(
    one_chain: Chain,
    original_sequence: &mut G,
    new_sequence: &mut G,
    original_vcf_writer: &mut W,
    new_vcf_writer: &mut W,
    sv_len: usize,
) -> Result<(), LiftOverError> {
    // Write VCF
    let mut reference_current = one_chain.original_start;
    let mut query_current = one_chain.new_start;
    for one_interval in one_chain.chain_interval {
        let reference_next = reference_current + one_interval.size;
        let query_next = query_current + one_interval.size;
        let (interval_reference_start, interval_reference_end) = convert_position(
            reference_current,
            reference_next,
            one_chain.original_chromosome.length,
            one_chain.original_strand,
        );
        let (interval_query_start, interval_query_end) = convert_position(
            query_current,
            query_next,
            one_chain.new_chromosome.length,
            one_chain.new_strand,
        );

        let reference_sequence_data = original_sequence.get_sequence(
            &one_chain.original_chromosome.name,
            interval_reference_start,
            interval_reference_end,
        )?;
        let query_sequence_data_tmp = new_sequence.get_sequence(
            &one_chain.new_chromosome.name,
            interval_query_start,
            interval_query_end,
        )?;
        let query_sequence_data = if one_chain.original_strand == one_chain.new_strand {
            query_sequence_data_tmp
        } else {
            reverse_complement(&query_sequence_data_tmp)
        };
        let different_positions: Vec<_> = reference_sequence_data
            .iter()
            .zip(query_sequence_data.iter())
            .enumerate()
            .filter(|(_, (x, y))| !compare_one_acid(**x, **y))
            .collect();

        // TODO: merge MNV
        assert_eq!(one_chain.original_strand, Strand::Forward);
        for one_position in different_positions {
            let query_vcf_pos = match one_chain.new_strand {
                Strand::Forward => one_position.0 as u64 + query_current + 1,
                Strand::Reverse => {
                    one_chain.new_chromosome.length - (one_position.0 as u64 + query_current + 1)
                        + 1
                }
            };

            write!(
                original_vcf_writer,
                "{}\t{}\t.\t",
                one_chain.original_chromosome.name,
                one_position.0 as u64 + reference_current + 1,
            )?;
            original_vcf_writer.write_all(&[*(one_position.1).0])?;
            original_vcf_writer.write_all(b"\t")?;
            original_vcf_writer.write_all(&[*(one_position.1).1])?;
            writeln!(
                original_vcf_writer,
                "\t.\t.\tTARGET_CHROM={};TARGET_POS={};CHAIN_ID={};STRAND={}",
                one_chain.new_chromosome.name,
                query_vcf_pos,
                one_chain.chain_id,
                one_chain.new_strand
            )?;

            // write query VCF
            write!(
                new_vcf_writer,
                "{}\t{}\t.\t",
                one_chain.new_chromosome.name, query_vcf_pos,
            )?;
            match one_chain.new_strand {
                Strand::Forward => {
                    new_vcf_writer.write_all(&[*(one_position.1).1])?;
                    new_vcf_writer.write_all(b"\t")?;
                    new_vcf_writer.write_all(&[*(one_position.1).0])?;
                }
                Strand::Reverse => {
                    new_vcf_writer.write_all(&[reverse_acid(*(one_position.1).1)])?;
                    new_vcf_writer.write_all(b"\t")?;
                    new_vcf_writer.write_all(&[reverse_acid(*(one_position.1).0)])?;
                }
            }
            writeln!(
                new_vcf_writer,
                "\t.\t.\tTARGET_CHROM={};TARGET_POS={};CHAIN_ID={};STRAND={}",
                one_chain.original_chromosome.name,
                one_position.0 as u64 + reference_current + 1,
                one_chain.chain_id,
                one_chain.new_strand
            )?;
        }

        reference_current = reference_next;
        query_current = query_next;

        // process gap
        let reference_next = reference_current + one_interval.difference_original.unwrap_or(0);
        let query_next = query_current + one_interval.difference_new.unwrap_or(0);
        if (reference_next == reference_current) || (query_next == query_current) {
            reference_current -= 1;
            query_current -= 1;
        }
        let (interval_reference_start, interval_reference_end) = convert_position(
            reference_current,
            reference_next,
            one_chain.original_chromosome.length,
            one_chain.original_strand,
        );
        let (interval_query_start, interval_query_end) = convert_position(
            query_current,
            query_next,
            one_chain.new_chromosome.length,
            one_chain.new_strand,
        );

        let reference_sequence_data = original_sequence.get_sequence(
            &one_chain.original_chromosome.name,
            interval_reference_start,
            interval_reference_end,
        )?;
        let query_sequence_data_tmp = new_sequence.get_sequence(
            &one_chain.new_chromosome.name,
            interval_query_start,
            interval_query_end,
        )?;
        let query_sequence_data = if one_chain.original_strand == one_chain.new_strand {
            query_sequence_data_tmp
        } else {
            reverse_complement(&query_sequence_data_tmp)
        };

        if reference_sequence_data != query_sequence_data {
            let query_vcf_pos = match one_chain.new_strand {
                Strand::Forward => query_current + 1,
                Strand::Reverse => {
                    one_chain.new_chromosome.length
                        - (query_current + query_sequence_data.len() as u64)
                        + 1
                }
            };
            if reference_sequence_data.len() > sv_len || query_sequence_data.len() > sv_len {
                // write as SV

                // write Reference VCF
                let (sv_type, sv_len) = if reference_sequence_data.len() == 1 {
                    ("INS", query_sequence_data.len())
                } else if query_sequence_data.len() == 1 {
                    ("DEL", reference_sequence_data.len())
                } else {
                    ("INDEL", reference_sequence_data.len())
                };

                write!(
                    original_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.original_chromosome.name,
                    reference_current + 1,
                )?;
                original_vcf_writer.write_all(&reference_sequence_data[0..1])?;
                writeln!(
                    original_vcf_writer,
                    "\t<{}>\t.\t.\tEND={};TARGET_CHROM={};TARGET_POS={};SVTYPE={};SVLEN={};CHAIN_ID={};STRAND={}",
                    sv_type,
                    reference_current + reference_sequence_data.len() as u64,
                    one_chain.new_chromosome.name,
                    query_vcf_pos,
                    sv_type,
                    sv_len,
                    one_chain.chain_id,
                    one_chain.new_strand
                )?;

                // Write Query VCF
                let sv_type_query = if reference_sequence_data.len() == 1 {
                    "DEL"
                } else if query_sequence_data.len() == 1 {
                    "INS"
                } else {
                    "INDEL"
                };

                write!(
                    new_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.new_chromosome.name, query_vcf_pos,
                )?;
                let query_vcf_ref = match one_chain.new_strand {
                    Strand::Forward => query_sequence_data,
                    Strand::Reverse => reverse_complement(&query_sequence_data),
                };
                new_vcf_writer.write_all(&query_vcf_ref[0..1])?;
                writeln!(
                    new_vcf_writer,
                    "\t<{}>\t.\t.\tEND={};TARGET_CHROM={};TARGET_POS={};SVTYPE={};SVLEN={};CHAIN_ID={};STRAND={}",
                    sv_type_query,
                    query_vcf_pos + query_vcf_ref.len() as u64,
                    one_chain.original_chromosome.name,
                    query_vcf_pos,
                    sv_type_query,
                    sv_len,
                    one_chain.chain_id,
                    one_chain.new_strand
                )?;
            } else {
                write!(
                    original_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.original_chromosome.name,
                    reference_current + 1,
                )?;
                original_vcf_writer.write_all(&reference_sequence_data)?;
                original_vcf_writer.write_all(b"\t")?;
                original_vcf_writer.write_all(&query_sequence_data)?;
                writeln!(
                    original_vcf_writer,
                    "\t.\t.\tTARGET_CHROM={};TARGET_POS={};CHAIN_ID={};STRAND={}",
                    one_chain.new_chromosome.name,
                    query_vcf_pos,
                    one_chain.chain_id,
                    one_chain.new_strand
                )?;

                // Write Query VCF
                write!(
                    new_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.new_chromosome.name, query_vcf_pos,
                )?;
                let query_vcf_ref = match one_chain.new_strand {
                    Strand::Forward => query_sequence_data,
                    Strand::Reverse => reverse_complement(&query_sequence_data),
                };
                new_vcf_writer.write_all(&query_vcf_ref)?;
                new_vcf_writer.write_all(b"\t")?;
                let query_vcf_alt = match one_chain.new_strand {
                    Strand::Forward => reference_sequence_data,
                    Strand::Reverse => reverse_complement(&reference_sequence_data),
                };
                new_vcf_writer.write_all(&query_vcf_alt)?;
                writeln!(
                    new_vcf_writer,
                    "\t.\t.\tTARGET_CHROM={};TARGET_POS={};CHAIN_ID={};STRAND={}",
                    one_chain.original_chromosome.name,
                    reference_current + 1,
                    one_chain.chain_id,
                    one_chain.new_strand
                )?;
            }
        }

        // TODO: write query VCF
        reference_current = reference_next;
        query_current = query_next;
    }
    Ok(())
}

fn convert_position(start: u64, stop: u64, length: u64, strand: Strand) -> (u64, u64) {
    match strand {
        Strand::Forward => (start, stop),
        Strand::Reverse => (length - stop, length - start),
    }
}

fn compare_one_acid(acid1: u8, acid2: u8) -> bool {
    let acid1_normed = match acid1 {
        b'a' | b'A' => b'A',
        b't' | b'T' => b'T',
        b'c' | b'C' => b'C',
        b'g' | b'G' => b'G',
        _ => acid1,
    };
    let acid2_normed = match acid2 {
        b'a' | b'A' => b'A',
        b't' | b'T' => b'T',
        b'c' | b'C' => b'C',
        b'g' | b'G' => b'G',
        _ => acid2,
    };
    acid1_normed == acid2_normed
}

fn write_vcf_header_for_sequence(
    writer: &mut dyn io::Write,
    sequence: &dyn GenomeSequence,
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

#[cfg(test)]
mod test {
    use clap::Parser;

    use crate::Cli;

    use super::*;
    use std::fs;

    #[test]
    fn test_chain_to_bed_vcf_helper_chr22() -> anyhow::Result<()> {
        fs::create_dir_all("../target/test-output/chain2bed")?;
        chain_to_bed_vcf_helper(
            "../liftover-rs/testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain",
            "../liftover-rs/testfiles/genomes/GRCh38/GRCh38.chr22.genome.fa",
            "../liftover-rs/testfiles/genomes/GRCh37/GRCh37.chr22.genome.fa",
            "../target/test-output/chain2bed/chain-to-bed-vcf--GRCh38.vcf",
            "../target/test-output/chain2bed/chain-to-bed-vcf--GRCh37.vcf",
            "../target/test-output/chain2bed/chain-to-bed-vcf--GRCh38.bed",
            "../target/test-output/chain2bed/chain-to-bed-vcf--GRCh37.bed",
            50,
        )?;

        // TODO: Check result

        Ok(())
    }

    #[test]
    fn test_chain_to_bed_vcf_helper_chr22_cli() -> anyhow::Result<()> {
        fs::create_dir_all("../target/test-output/chain2bed-2")?;

        let cli = Cli::parse_from(&[
            "transanno",
            "chain-to-bed-vcf",
            "../liftover-rs/testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain",
            "--original",
            "../liftover-rs/testfiles/genomes/GRCh38/GRCh38.chr22.genome.fa",
            "--new",
            "../liftover-rs/testfiles/genomes/GRCh37/GRCh37.chr22.genome.fa",
            "--output-new-bed",
            "../target/test-output/chain2bed-2/new-GRCh37.bed",
            "--output-new-vcf",
            "../target/test-output/chain2bed-2/new-GRCh37.vcf.gz",
            "--output-original-bed",
            "../target/test-output/chain2bed-2/original-GRCh38.bed",
            "--output-original-vcf",
            "../target/test-output/chain2bed-2/original-GRCh38.vcf.gz",
        ]);
        cli.command.run()?;

        // TODO: Check result

        Ok(())
    }
}
