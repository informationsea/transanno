use super::*;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct JsonWrapper {
    list: Vec<ExpectedRegion>,
}

#[derive(Debug, PartialEq, Eq, Hash, Deserialize, Clone)]
struct ExpectedRegion {
    chrom: String,
    start: u64,
    end: u64,
    mapped: Vec<MappedRegion>,
}

#[derive(Debug, PartialEq, Eq, Hash, Deserialize, Clone)]
struct MappedRegion {
    chrom: String,
    start: u64,
    end: u64,
    strand: String,
}

#[test]
fn test_ucsc_result() -> anyhow::Result<()> {
    let lift_over = PositionLiftOver::load(autocompress::open(
        "./testfiles/ucsc/hg19ToHg38.over.chain.gz",
    )?)?;

    let expected: JsonWrapper = serde_json::from_reader(autocompress::open(
        "./testfiles/ucsc/hg19-regions-mapped-to-hg38.json",
    )?)?;

    Ok(())
}
