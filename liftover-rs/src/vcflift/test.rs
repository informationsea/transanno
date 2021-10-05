use super::*;
use crate::chain::ChainFile;
use bio::io::fasta::IndexedReader;
use std::fs::File;
use std::str;

fn new_vcf_liftover(
    parameter: VCFLiftOverParameters,
) -> anyhow::Result<VCFLiftOver<IndexedReader<File>>> {
    let chain = ChainFile::load(
        &include_bytes!("../../testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain")[..],
    )?;
    let mut grch37 = IndexedReader::from_file(&"testfiles/genomes/GRCh37/GRCh37.chr22.genome.fa")?;
    let mut grch38 = IndexedReader::from_file(&"testfiles/genomes/GRCh38/GRCh38.chr22.genome.fa")?;
    Ok(VCFLiftOver::new(
        VariantLiftOver::new(chain.left_align(&mut grch38, &mut grch37)?, grch38, grch37),
        parameter,
    ))
}

#[test]
fn test_header_lift() -> anyhow::Result<()> {
    let mut sample_header = &include_bytes!("testfiles/test-header.vcf")[..];
    let reader = VCFReader::new(&mut sample_header)?;
    let vcf_lift = new_vcf_liftover(VCFLiftOverParameters::new())?;

    let (_lifted_header, rewrite_target) = vcf_lift.lift_header(&reader.header)?;

    let expected = VCFHeaderRewriteTarget {
        info_alt: vec![b"number_a".to_vec()].into_iter().collect(),
        info_ref: vec![b"number_r".to_vec()].into_iter().collect(),
        info_genotype: vec![b"number_g".to_vec()].into_iter().collect(),
        format_alt: vec![b"f_number_a".to_vec()].into_iter().collect(),
        format_ref: vec![b"f_number_r".to_vec()].into_iter().collect(),
        format_genotype: vec![b"f_number_g".to_vec()].into_iter().collect(),
        allele_count: vec![
            b"AC".to_vec(),
            b"AC_nfe_seu".to_vec(),
            b"non_topmed_AC_amr".to_vec(),
            b"non_neuro_AC".to_vec(),
        ]
        .into_iter()
        .collect(),
        allele_number: vec![
            b"AN".to_vec(),
            b"AN_nfe_seu".to_vec(),
            b"non_topmed_AN_amr".to_vec(),
            b"non_neuro_AN".to_vec(),
        ]
        .into_iter()
        .collect(),
        allele_count_to_allele_number: vec![
            (b"AC".to_vec(), b"AN".to_vec()),
            (b"AC_nfe_seu".to_vec(), b"AN_nfe_seu".to_vec()),
            (b"non_topmed_AC_amr".to_vec(), b"non_topmed_AN_amr".to_vec()),
            (b"non_neuro_AC".to_vec(), b"non_neuro_AN".to_vec()),
        ]
        .into_iter()
        .collect(),
        allele_frequency: vec![
            b"AF".to_vec(),
            b"AF_nfe_seu".to_vec(),
            b"non_topmed_AF_amr".to_vec(),
            b"non_neuro_AF".to_vec(),
        ]
        .into_iter()
        .collect(),
        format_gt: false,
    };
    assert_eq!(rewrite_target, expected);

    Ok(())
}

#[test]
fn test_lift_noswap_vcf() -> anyhow::Result<()> {
    let mut reader = VCFReader::new(&include_bytes!("testfiles/original.vcf")[..])?;
    let mut expected = VCFReader::new(&include_bytes!("testfiles/mapped-noswap.vcf")[..])?;
    let mut vcf_lift = new_vcf_liftover(VCFLiftOverParameters::new().do_not_swap_ref_alt(true))?;

    let (_lifted_header, rewrite_target) = vcf_lift.lift_header(&reader.header)?;

    while let Some(record) = reader.next_record()? {
        let lifted_record = vcf_lift.lift_record(&record, &rewrite_target)?;
        let expected_record = expected.next_record()?.unwrap();
        match lifted_record {
            VCFLiftOverResult::Succeeded(succeeded_records) => {
                assert_eq!(succeeded_records.len(), 1);
                let mut expected_bytes: Vec<u8> = Vec::new();
                expected_record.write(&mut expected_bytes)?;
                let mut lifted_bytes: Vec<u8> = Vec::new();
                succeeded_records[0].write(&mut lifted_bytes)?;
                assert_eq!(
                    str::from_utf8(&expected_bytes).unwrap(),
                    str::from_utf8(&lifted_bytes).unwrap()
                );
            }
            VCFLiftOverResult::Failed(_) => panic!(),
        }
    }

    Ok(())
}

#[test]
fn test_lift_swap_vcf() -> anyhow::Result<()> {
    let mut reader = VCFReader::new(&include_bytes!("testfiles/original.vcf")[..])?;
    let mut expected = VCFReader::new(&include_bytes!("testfiles/mapped-swap.vcf")[..])?;
    let mut vcf_lift = new_vcf_liftover(VCFLiftOverParameters::new())?;

    let (_lifted_header, rewrite_target) = vcf_lift.lift_header(&reader.header)?;

    while let Some(record) = reader.next_record()? {
        let lifted_record = vcf_lift.lift_record(&record, &rewrite_target)?;
        let expected_record = expected.next_record()?.unwrap();
        match lifted_record {
            VCFLiftOverResult::Succeeded(succeeded_records) => {
                assert_eq!(succeeded_records.len(), 1);
                let mut expected_bytes: Vec<u8> = Vec::new();
                expected_record.write(&mut expected_bytes)?;
                let mut lifted_bytes: Vec<u8> = Vec::new();
                succeeded_records[0].write(&mut lifted_bytes)?;
                assert_eq!(
                    str::from_utf8(&expected_bytes).unwrap(),
                    str::from_utf8(&lifted_bytes).unwrap()
                );
            }
            VCFLiftOverResult::Failed(_) => panic!(),
        }
    }

    Ok(())
}
