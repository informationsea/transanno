use std::{fs::File, io::BufWriter};

use anyhow::Context;
use autocompress::autodetect_open;
use bio::io::fasta::IndexedReader;

use super::*;

fn try_left_align(
    chain_path: &str,
    output_path: &str,
    original_fasta_path: &str,
    new_fasta_path: &str,
) -> anyhow::Result<()> {
    let mut original_fasta =
        IndexedReader::from_file(&original_fasta_path).context("Failed to load original fasta")?;
    let mut new_fasta =
        IndexedReader::from_file(&new_fasta_path).context("Failed to load new fasta")?;
    let chain = ChainFile::load(autodetect_open(chain_path).context("Failed to open chain file")?)
        .context("Failed to load chain file")?;
    eprintln!("chain file loaded: {}", chain.chain_list.len());
    chain
        .left_align(&mut original_fasta, &mut new_fasta)
        .context("Failed to run left align")?;
    let mut writer =
        BufWriter::new(File::create(output_path).context("Failed to open chain file writer")?);
    chain
        .write(&mut writer)
        .context("Failed to write chain file")?;
    Ok(())
}

#[test]
fn left_align_hg19_to_hg38() -> anyhow::Result<()> {
    try_left_align(
        "testfiles/ucsc/hg19ToHg38.over.chain.gz",
        "../target/chain-ucsc-hg19ToHg38-left-align.chain",
        "testfiles/ucsc/hg19.fa",
        "testfiles/ucsc/hg38.fa",
    )?;
    Ok(())
}

#[test]
fn left_align_hg38_to_hg19() -> anyhow::Result<()> {
    try_left_align(
        "testfiles/ucsc/hg38ToHg19.over.chain.gz",
        "../target/chain-ucsc-hg38ToHg19-left-align.chain",
        "testfiles/ucsc/hg38.fa",
        "testfiles/ucsc/hg19.fa",
    )?;
    Ok(())
}

#[test]
fn left_align_hg38_to_hg19_minimap2() -> anyhow::Result<()> {
    try_left_align(
        "testfiles/minimap2/hg38.p14__to__hg19.p13.plusMT.chain.gz",
        "../target/chain-ucsc-hg38ToHg19-left-align-minimap2.chain.gz",
        "testfiles/ucsc/hg38.fa",
        "testfiles/ucsc/hg19.fa",
    )?;
    Ok(())
}

#[test]
fn left_align_hg38_to_hs1() -> anyhow::Result<()> {
    try_left_align(
        "testfiles/ucsc/hg38ToHs1.over.chain.gz",
        "../target/chain-ucsc-hg38ToHs1-left-align.chain",
        "testfiles/ucsc/hg38.fa",
        "testfiles/ucsc/hs1.fa",
    )?;
    Ok(())
}

#[test]
fn left_align_hg38_to_hs1_minimap2() -> anyhow::Result<()> {
    try_left_align(
        "testfiles/minimap2/hg38.p14__to__hs1.chain.gz",
        "../target/chain-ucsc-hg38ToHs1-left-align.chain",
        "testfiles/ucsc/hg38.fa",
        "testfiles/ucsc/hs1.fa",
    )?;
    Ok(())
}

#[test]
fn left_align_hg38_to_mm39() -> anyhow::Result<()> {
    try_left_align(
        "testfiles/ucsc/hg38ToMm39.over.chain.gz",
        "../target/chain-ucsc-hg38ToMm39-left-align.chain",
        "testfiles/ucsc/hg38.fa",
        "testfiles/ucsc/mm39.fa",
    )?;
    Ok(())
}
