use crate::chain::{Chain, Strand};
use crate::defs::{
    adaptive_create, adaptive_open, reverse_acid, reverse_complement, GenomeSequence,
};
use crate::LiftOverError;
use bio::io::fasta::IndexedReader;
use clap::ArgMatches;
use std::io;

pub fn chain_to_bed_vcf(matches: &ArgMatches) {
    chain_to_bed_vcf_helper(
        matches.value_of("chain").unwrap(),
        matches.value_of("reference-sequence").unwrap(),
        matches.value_of("query-sequence").unwrap(),
        matches.value_of("reference-vcf").unwrap(),
        matches.value_of("query-vcf").unwrap(),
        matches.value_of("reference-bed").unwrap(),
        matches.value_of("query-bed").unwrap(),
        matches.value_of("svlen").unwrap().parse().unwrap(),
    )
    .expect("failed to create BED and VCF");
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
        crate::chain::ChainFile::load(adaptive_open(chain_path).expect("Cannot open chain file"))
            .expect("Cannot parse chain file");

    let mut query_vcf_writer =
        adaptive_create(query_vcf_path).expect("Cannot create Query VCF file");
    write_vcf_header_for_sequence(&mut query_vcf_writer, &query_sequence, query_sequence_path)
        .expect("Cannot write query VCF");
    let mut query_bed_writer =
        adaptive_create(query_bed_path).expect("Cannot create Query BED file");

    let mut reference_vcf_writer =
        adaptive_create(reference_vcf_path).expect("Cannot create reference VCF file");
    write_vcf_header_for_sequence(
        &mut reference_vcf_writer,
        &reference_sequence,
        reference_sequence_path,
    )
    .expect("Cannot write reference VCF");

    let mut reference_bed_writer =
        adaptive_create(reference_bed_path).expect("Cannot create reference BED file");

    for one_chain in chain_file.chain_list {
        // write BED
        let (reference_start, reference_end) = convert_position(
            one_chain.reference_start,
            one_chain.reference_end,
            one_chain.reference_chromosome.length,
            one_chain.reference_strand,
        );
        let (query_start, query_end) = convert_position(
            one_chain.query_start,
            one_chain.query_end,
            one_chain.query_chromosome.length,
            one_chain.query_strand,
        );

        writeln!(
            reference_bed_writer,
            "{}\t{}\t{}\tchain_id:{};{}:{}-{}\t{}\t{}",
            one_chain.reference_chromosome.name,
            reference_start,
            reference_end,
            one_chain.chain_id,
            one_chain.query_chromosome.name,
            query_start + 1,
            query_end,
            one_chain.score,
            if one_chain.reference_strand == one_chain.query_strand {
                "+"
            } else {
                "-"
            },
        )
        .expect("Cannot write reference BED");

        writeln!(
            query_bed_writer,
            "{}\t{}\t{}\tchain_id:{};{}:{}-{}\t{}\t{}",
            one_chain.query_chromosome.name,
            query_start,
            query_end,
            one_chain.chain_id,
            one_chain.reference_chromosome.name,
            reference_start + 1,
            reference_end,
            one_chain.score,
            if one_chain.reference_strand == one_chain.query_strand {
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
    let mut reference_current = one_chain.reference_start;
    let mut query_current = one_chain.query_start;
    for one_interval in one_chain.chain_interval {
        let reference_next = reference_current + one_interval.size;
        let query_next = query_current + one_interval.size;
        let (interval_reference_start, interval_reference_end) = convert_position(
            reference_current,
            reference_next,
            one_chain.reference_chromosome.length,
            one_chain.reference_strand,
        );
        let (interval_query_start, interval_query_end) = convert_position(
            query_current,
            query_next,
            one_chain.query_chromosome.length,
            one_chain.query_strand,
        );

        let reference_sequence_data = reference_sequence.get_sequence(
            &one_chain.reference_chromosome.name,
            interval_reference_start,
            interval_reference_end,
        )?;
        let query_sequence_data_tmp = query_sequence.get_sequence(
            &one_chain.query_chromosome.name,
            interval_query_start,
            interval_query_end,
        )?;
        let query_sequence_data = if one_chain.reference_strand == one_chain.query_strand {
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
        assert_eq!(one_chain.reference_strand, Strand::Forward);
        for one_position in different_positions {
            let query_vcf_pos = match one_chain.query_strand {
                Strand::Forward => one_position.0 as u64 + query_current + 1,
                Strand::Reverse => {
                    one_chain.query_chromosome.length - (one_position.0 as u64 + query_current + 1)
                        + 1
                }
            };

            write!(
                reference_vcf_writer,
                "{}\t{}\t.\t",
                one_chain.reference_chromosome.name,
                one_position.0 as u64 + reference_current + 1,
            )?;
            reference_vcf_writer.write_all(&[*(one_position.1).0])?;
            reference_vcf_writer.write_all(b"\t")?;
            reference_vcf_writer.write_all(&[*(one_position.1).1])?;
            writeln!(
                reference_vcf_writer,
                "\t.\t.\tTARGET_CHROM={};TARGET_POS={}",
                one_chain.query_chromosome.name, query_vcf_pos
            )?;

            // write query VCF
            write!(
                query_vcf_writer,
                "{}\t{}\t.\t",
                one_chain.query_chromosome.name, query_vcf_pos,
            )?;
            match one_chain.query_strand {
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
                "\t.\t.\tTARGET_CHROM={};TARGET_POS={}",
                one_chain.reference_chromosome.name,
                one_position.0 as u64 + reference_current + 1
            )?;
        }

        reference_current = reference_next;
        query_current = query_next;

        // process gap
        let reference_next = reference_current + one_interval.difference_reference.unwrap_or(0);
        let query_next = query_current + one_interval.difference_query.unwrap_or(0);
        if (reference_next == reference_current) || (query_next == query_current) {
            reference_current -= 1;
            query_current -= 1;
        }
        let (interval_reference_start, interval_reference_end) = convert_position(
            reference_current,
            reference_next,
            one_chain.reference_chromosome.length,
            one_chain.reference_strand,
        );
        let (interval_query_start, interval_query_end) = convert_position(
            query_current,
            query_next,
            one_chain.query_chromosome.length,
            one_chain.query_strand,
        );

        let reference_sequence_data = reference_sequence.get_sequence(
            &one_chain.reference_chromosome.name,
            interval_reference_start,
            interval_reference_end,
        )?;
        let query_sequence_data_tmp = query_sequence.get_sequence(
            &one_chain.query_chromosome.name,
            interval_query_start,
            interval_query_end,
        )?;
        let query_sequence_data = if one_chain.reference_strand == one_chain.query_strand {
            query_sequence_data_tmp
        } else {
            reverse_complement(&query_sequence_data_tmp)
        };

        if reference_sequence_data != query_sequence_data {
            let query_vcf_pos = match one_chain.query_strand {
                Strand::Forward => query_current + 1,
                Strand::Reverse => {
                    one_chain.query_chromosome.length
                        - (query_current + query_sequence_data.len() as u64)
                        + 1
                }
            };
            if reference_sequence_data.len() > sv_len || query_sequence_data.len() > sv_len {
                // write as SV

                // write Reference VCF
                let sv_type = if reference_sequence_data.len() == 1 {
                    "INS"
                } else if query_sequence_data.len() == 1 {
                    "DEL"
                } else {
                    "INDEL"
                };

                write!(
                    reference_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.reference_chromosome.name,
                    reference_current + 1,
                )?;
                reference_vcf_writer.write_all(&reference_sequence_data[0..1])?;
                writeln!(
                    reference_vcf_writer,
                    "\t<{}>\t.\t.\tEND={};TARGET_CHROM={};TARGET_POS={};SVTYPE={};SVLEN={}",
                    sv_type,
                    reference_current + reference_sequence_data.len() as u64,
                    one_chain.query_chromosome.name,
                    query_vcf_pos,
                    sv_type,
                    query_sequence_data.len()
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
                    one_chain.query_chromosome.name, query_vcf_pos,
                )?;
                let query_vcf_ref = match one_chain.query_strand {
                    Strand::Forward => query_sequence_data,
                    Strand::Reverse => reverse_complement(&query_sequence_data),
                };
                query_vcf_writer.write_all(&query_vcf_ref[0..1])?;
                writeln!(
                    query_vcf_writer,
                    "\t<{}>\t.\t.\tEND={};TARGET_CHROM={};TARGET_POS={};SVTYPE={};SVLEN={}",
                    sv_type_query,
                    query_current + query_vcf_ref.len() as u64,
                    one_chain.reference_chromosome.name,
                    query_vcf_pos,
                    sv_type_query,
                    reference_sequence_data.len()
                )?;
            } else {
                write!(
                    reference_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.reference_chromosome.name,
                    reference_current + 1,
                )?;
                reference_vcf_writer.write_all(&reference_sequence_data)?;
                reference_vcf_writer.write_all(b"\t")?;
                reference_vcf_writer.write_all(&query_sequence_data)?;
                writeln!(
                    reference_vcf_writer,
                    "\t.\t.\tTARGET_CHROM={};TARGET_POS={}",
                    one_chain.query_chromosome.name, query_vcf_pos
                )?;

                // Write Query VCF
                write!(
                    query_vcf_writer,
                    "{}\t{}\t.\t",
                    one_chain.query_chromosome.name, query_vcf_pos,
                )?;
                let query_vcf_ref = match one_chain.query_strand {
                    Strand::Forward => query_sequence_data,
                    Strand::Reverse => reverse_complement(&query_sequence_data),
                };
                query_vcf_writer.write_all(&query_vcf_ref)?;
                query_vcf_writer.write_all(b"\t")?;
                let query_vcf_alt = match one_chain.query_strand {
                    Strand::Forward => reference_sequence_data,
                    Strand::Reverse => reverse_complement(&reference_sequence_data),
                };
                query_vcf_writer.write_all(&query_vcf_alt)?;
                writeln!(
                    query_vcf_writer,
                    "\t.\t.\tTARGET_CHROM={};TARGET_POS={}",
                    one_chain.reference_chromosome.name,
                    reference_current + 1
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
    fn test_chain_to_bed_vcf_helper_test_seq() -> Result<(), LiftOverError> {
        fs::create_dir_all("target/test-output/chain2bed")?;

        chain_to_bed_vcf_helper(
            "testfiles/lift-variant-test/sequence/seq-a--to--seq-b.chain.gz",
            "testfiles/lift-variant-test/sequence/seq-a.fa",
            "testfiles/lift-variant-test/sequence/seq-b.fa",
            "target/test-output/chain2bed/chain-to-bed-vcf--seq-a.vcf",
            "target/test-output/chain2bed/chain-to-bed-vcf--seq-b.vcf",
            "target/test-output/chain2bed/chain-to-bed-vcf--seq-a.bed",
            "target/test-output/chain2bed/chain-to-bed-vcf--seq-b.bed",
            50,
        )?;

        // TODO: Check result

        Ok(())
    }

    #[test]
    fn test_chain_to_bed_vcf_helper_chr22() -> Result<(), LiftOverError> {
        fs::create_dir_all("target/test-output/chain2bed")?;
        chain_to_bed_vcf_helper(
            "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
            "testfiles/genome/hg19-chr22.fa",
            "testfiles/genome/hg38-chr22.fa",
            "target/test-output/chain2bed/chain-to-bed-vcf--hg19.vcf",
            "target/test-output/chain2bed/chain-to-bed-vcf--hg38.vcf",
            "target/test-output/chain2bed/chain-to-bed-vcf--hg19.bed",
            "target/test-output/chain2bed/chain-to-bed-vcf--hg38.bed",
            50,
        )?;

        // TODO: Check result

        Ok(())
    }
}
