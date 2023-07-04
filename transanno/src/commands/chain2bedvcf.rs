use autocompress::{create, open, CompressionLevel};
use bio::io::fasta::IndexedReader;
use clap::Args;
use liftover::chain::{Chain, Strand};
use liftover::LiftOverError;
use liftover::{reverse_acid, reverse_complement, GenomeSequence};
use log::warn;
use std::io::{self, Write};

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
    reference_bed: String,
    #[arg(
        long = "output-new-bed",
        alias = "output-query-bed",
        short = 'd',
        help = "Output new assembly BED file (Not sorted)"
    )]
    query_bed: String,
    #[arg(
        long = "output-original-vcf",
        alias = "output-reference-vcf",
        short = 'v',
        help = "Output original assembly VCF file (Not sorted)"
    )]
    reference_vcf: String,
    #[arg(
        long = "output-new-vcf",
        alias = "output-query-vcf",
        short = 'c',
        help = "Output new assembly VCF file (Not sorted)"
    )]
    query_vcf: String,
    #[arg(
        long = "original",
        alias = "reference",
        short = 'r',
        help = "Original assembly FASTA (.fai file is required)"
    )]
    reference_sequence: String,
    #[arg(
        long = "new",
        alias = "query",
        short = 'q',
        help = "New assembly FASTA (.fai file is required)"
    )]
    query_sequence: String,
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
            &self.reference_sequence,
            &self.query_sequence,
            &self.reference_vcf,
            &self.query_vcf,
            &self.reference_bed,
            &self.query_bed,
            self.svlen,
        )?;
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn chain_to_bed_vcf_helper(
    chain_path: &str,
    reference_sequence_path: &str,
    query_sequence_path: &str,
    reference_vcf_path: &str,
    query_vcf_path: &str,
    reference_bed_path: &str,
    query_bed_path: &str,
    sv_len: usize,
) -> Result<(), LiftOverError> {
    let mut reference_sequence =
        IndexedReader::from_file(&reference_sequence_path).expect("Cannot open reference sequence");
    let mut query_sequence =
        IndexedReader::from_file(&query_sequence_path).expect("Cannot open query sequence");
    let chain_file =
        liftover::chain::ChainFile::load(open(chain_path).expect("Cannot open chain file"))
            .expect("Cannot parse chain file");

    let mut query_vcf_writer =
        create(query_vcf_path, CompressionLevel::Default).expect("Cannot create Query VCF file");
    write_vcf_header_for_sequence(&mut query_vcf_writer, &query_sequence, query_sequence_path)
        .expect("Cannot write query VCF");
    let mut query_bed_writer =
        create(query_bed_path, CompressionLevel::Default).expect("Cannot create Query BED file");

    let mut reference_vcf_writer = create(reference_vcf_path, CompressionLevel::Default)
        .expect("Cannot create reference VCF file");
    write_vcf_header_for_sequence(
        &mut reference_vcf_writer,
        &reference_sequence,
        reference_sequence_path,
    )
    .expect("Cannot write reference VCF");

    let mut reference_bed_writer = create(reference_bed_path, CompressionLevel::Default)
        .expect("Cannot create reference BED file");

    for one_chain in chain_file.chain_list {
        // Check chromosome length
        match one_chain.check_sequence_consistency(&mut reference_sequence, &mut query_sequence) {
            Ok(_) => (),
            Err(LiftOverError::ChromosomeNotFound(_)) => {
                warn!("Skip chain ID: {}", one_chain.chain_id);
                continue;
            }
            Err(e) => return Err(e),
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
            reference_bed_writer,
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
            query_bed_writer,
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
            &mut reference_sequence,
            &mut query_sequence,
            &mut reference_vcf_writer,
            &mut query_vcf_writer,
            sv_len,
        )?;
    }

    Ok(())
}

fn write_vcf_entry<G: GenomeSequence, W: io::Write>(
    one_chain: Chain,
    reference_sequence: &mut G,
    query_sequence: &mut G,
    reference_vcf_writer: &mut W,
    query_vcf_writer: &mut W,
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

        let reference_sequence_data = reference_sequence.get_sequence(
            &one_chain.original_chromosome.name,
            interval_reference_start,
            interval_reference_end,
        )?;
        let query_sequence_data_tmp = query_sequence.get_sequence(
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
                reference_vcf_writer,
                "{}\t{}\t.\t",
                one_chain.original_chromosome.name,
                one_position.0 as u64 + reference_current + 1,
            )?;
            reference_vcf_writer.write_all(&[*(one_position.1).0])?;
            reference_vcf_writer.write_all(b"\t")?;
            reference_vcf_writer.write_all(&[*(one_position.1).1])?;
            writeln!(
                reference_vcf_writer,
                "\t.\t.\tTARGET_CHROM={};TARGET_POS={};CHAIN_ID={};STRAND={}",
                one_chain.new_chromosome.name,
                query_vcf_pos,
                one_chain.chain_id,
                one_chain.new_strand
            )?;

            // write query VCF
            write!(
                query_vcf_writer,
                "{}\t{}\t.\t",
                one_chain.new_chromosome.name, query_vcf_pos,
            )?;
            match one_chain.new_strand {
                Strand::Forward => {
                    query_vcf_writer.write_all(&[*(one_position.1).1])?;
                    query_vcf_writer.write_all(b"\t")?;
                    query_vcf_writer.write_all(&[*(one_position.1).0])?;
                }
                Strand::Reverse => {
                    query_vcf_writer.write_all(&[reverse_acid(*(one_position.1).1)])?;
                    query_vcf_writer.write_all(b"\t")?;
                    query_vcf_writer.write_all(&[reverse_acid(*(one_position.1).0)])?;
                }
            }
            writeln!(
                query_vcf_writer,
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

        let reference_sequence_data = reference_sequence.get_sequence(
            &one_chain.original_chromosome.name,
            interval_reference_start,
            interval_reference_end,
        )?;
        let query_sequence_data_tmp = query_sequence.get_sequence(
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
                    reference_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.original_chromosome.name,
                    reference_current + 1,
                )?;
                reference_vcf_writer.write_all(&reference_sequence_data[0..1])?;
                writeln!(
                    reference_vcf_writer,
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
                    query_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.new_chromosome.name, query_vcf_pos,
                )?;
                let query_vcf_ref = match one_chain.new_strand {
                    Strand::Forward => query_sequence_data,
                    Strand::Reverse => reverse_complement(&query_sequence_data),
                };
                query_vcf_writer.write_all(&query_vcf_ref[0..1])?;
                writeln!(
                    query_vcf_writer,
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
                    reference_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.original_chromosome.name,
                    reference_current + 1,
                )?;
                reference_vcf_writer.write_all(&reference_sequence_data)?;
                reference_vcf_writer.write_all(b"\t")?;
                reference_vcf_writer.write_all(&query_sequence_data)?;
                writeln!(
                    reference_vcf_writer,
                    "\t.\t.\tTARGET_CHROM={};TARGET_POS={};CHAIN_ID={};STRAND={}",
                    one_chain.new_chromosome.name,
                    query_vcf_pos,
                    one_chain.chain_id,
                    one_chain.new_strand
                )?;

                // Write Query VCF
                write!(
                    query_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.new_chromosome.name, query_vcf_pos,
                )?;
                let query_vcf_ref = match one_chain.new_strand {
                    Strand::Forward => query_sequence_data,
                    Strand::Reverse => reverse_complement(&query_sequence_data),
                };
                query_vcf_writer.write_all(&query_vcf_ref)?;
                query_vcf_writer.write_all(b"\t")?;
                let query_vcf_alt = match one_chain.new_strand {
                    Strand::Forward => reference_sequence_data,
                    Strand::Reverse => reverse_complement(&reference_sequence_data),
                };
                query_vcf_writer.write_all(&query_vcf_alt)?;
                writeln!(
                    query_vcf_writer,
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
    use super::*;
    use std::fs;

    #[test]
    fn test_chain_to_bed_vcf_helper_chr22() -> Result<(), LiftOverError> {
        fs::create_dir_all("../target/test-output/chain2bed")?;
        chain_to_bed_vcf_helper(
            "../liftover-rs/testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain",
            "../liftover-rs/testfiles/genomes/GRCh38/GRCh38.chr22.genome.fa",
            "../liftover-rs/testfiles/genomes/GRCh37/GRCh37.chr22.genome.fa",
            "../target/test-output/chain2bed/chain-to-bed-vcf--GRCh37.vcf",
            "../target/test-output/chain2bed/chain-to-bed-vcf--GRCh38.vcf",
            "../target/test-output/chain2bed/chain-to-bed-vcf--GRCh37.bed",
            "../target/test-output/chain2bed/chain-to-bed-vcf--GRCh38.bed",
            50,
        )?;

        // TODO: Check result

        Ok(())
    }
}
