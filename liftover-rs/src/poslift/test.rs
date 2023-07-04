use super::*;

use crate::LiftOverError;

#[test]
fn test_target_region() {
    let value = TargetRegion {
        original_chromosome_index: 1,
        new_chromosome_index: 1,
        chain_index: 1,
        original_start: 10,
        original_end: 20,
        new_start: 30,
        new_end: 40,
        strand: Strand::Forward,
        is_in_gap: true,
    };
    assert_eq!(value.contains_original_position(9), false);
    assert_eq!(value.contains_original_position(10), true);
    assert_eq!(value.contains_original_position(19), true);
    assert_eq!(value.contains_original_position(20), false);

    assert_eq!(value.contains_new_position(29), false);
    assert_eq!(value.contains_new_position(30), true);
    assert_eq!(value.contains_new_position(39), true);
    assert_eq!(value.contains_new_position(40), false);

    assert_eq!(value.is_indel(), false);
}

#[test]
fn test_target_region2() {
    let value = TargetRegion {
        original_chromosome_index: 1,
        new_chromosome_index: 1,
        chain_index: 1,
        original_start: 10,
        original_end: 20,
        new_start: 30,
        new_end: 50,
        strand: Strand::Forward,
        is_in_gap: true,
    };
    assert_eq!(value.contains_original_position(9), false);
    assert_eq!(value.contains_original_position(10), true);
    assert_eq!(value.contains_original_position(19), true);
    assert_eq!(value.contains_original_position(20), false);

    assert_eq!(value.contains_new_position(29), false);
    assert_eq!(value.contains_new_position(30), true);
    assert_eq!(value.contains_new_position(40), true);
    assert_eq!(value.contains_new_position(39), true);
    assert_eq!(value.contains_new_position(50), false);

    assert_eq!(value.is_indel(), true);
}

use serde::Deserialize;
use std::collections::HashSet;

#[derive(Debug, Deserialize)]
struct JsonWrapper {
    list: Vec<ExpectedLiftOver>,
}

#[derive(Debug, PartialEq, Deserialize)]
struct ExpectedLiftOver {
    original_chrom: String,
    original_pos: u64,
    mapped: Vec<SimpleResult>,
}

#[derive(Debug, PartialEq, Deserialize, Eq, Hash, Clone)]
struct SimpleResult {
    chrom: String,
    pos: u64,
    strand: String,
}

#[allow(clippy::unreadable_literal)]
#[allow(clippy::cognitive_complexity)]
#[test]
fn test_liftover1() -> Result<(), LiftOverError> {
    let lift_over = PositionLiftOver::load(
        &include_bytes!("../../testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain")[..],
    )?;
    let expected_loader: JsonWrapper = serde_json::from_reader(
        &include_bytes!("../../testfiles/genomes/chain/GRCh38-to-GRCh37.test-position.json")[..],
    )?;

    for one in expected_loader.list {
        let mapped: HashSet<_> = lift_over
            .lift_position(&one.original_chrom, one.original_pos)
            .iter()
            .map(|x| SimpleResult {
                chrom: x.chromosome.name.to_string(),
                pos: x.position,
                strand: x.strand.to_string(),
            })
            .collect();
        let expected: HashSet<_> = one.mapped.iter().cloned().collect();
        assert_eq!(mapped, expected);
    }

    Ok(())
}

#[derive(Debug, PartialEq, Deserialize, Clone)]
struct ExpectedLiftRegionResult {
    lift: Vec<ExpectedRegion>,
}

#[derive(Debug, PartialEq, Eq, Hash, Deserialize, Clone)]
struct ExpectedRegion {
    original_chrom: String,
    original_start: u64,
    original_end: u64,
    original_strand: String,
    mapped: Vec<MappedRegion>,
}

#[derive(Debug, PartialEq, Eq, Hash, Deserialize, Clone)]
struct MappedRegion {
    chrom: String,
    start: u64,
    end: u64,
    strand: String,
}

#[allow(clippy::unreadable_literal)]
#[test]
fn test_lift_region() -> Result<(), LiftOverError> {
    let lift_over = PositionLiftOver::load(
        &include_bytes!("../../testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain")[..],
    )?;
    let expected_regions: ExpectedLiftRegionResult = serde_json::from_reader(
        &include_bytes!("../../testfiles/genomes/chain/mapped-region/GRCh38-to-GRCh37.test-region.GRCh38.mapped.GRCh37.json")[..],
    )?;

    // these regions contains large deletion.
    // UCSC LiftOver tool cannot map them.
    let start_blacklist: HashSet<_> = [
        11142000, 16164000, 17409000, 22827000, 24000000, 24063000, 49383000, 49386000,
    ]
    .iter()
    .collect();

    for expected in expected_regions.lift {
        if start_blacklist.contains(&expected.original_start) {
            continue;
        }
        let result = lift_over.lift_region(
            &expected.original_chrom,
            expected.original_start..expected.original_end,
        );
        let position_mapped: HashSet<_> = lift_over
            .lift_position(&expected.original_chrom, expected.original_start)
            .iter()
            .map(|x| match x.strand {
                Strand::Forward => x.position,
                Strand::Reverse => x.position + 1,
            })
            .chain(
                lift_over
                    .lift_position(&expected.original_chrom, expected.original_end - 1)
                    .iter()
                    .map(|x| match x.strand {
                        Strand::Forward => x.position + 1,
                        Strand::Reverse => x.position,
                    }),
            )
            .collect();
        let result_mapped: HashSet<_> = result
            .iter()
            .map(|x| MappedRegion {
                chrom: x.chromosome.name.to_string(),
                start: x.start,
                end: x.end,
                strand: match x.strand {
                    Strand::Forward => "+".to_string(),
                    Strand::Reverse => "-".to_string(),
                },
            })
            .collect();

        // start position and end position should be included in a chain
        let expected_mapped: HashSet<_> = expected
            .mapped
            .iter()
            .filter(|x| position_mapped.contains(&x.start) && position_mapped.contains(&x.end))
            .cloned()
            .collect();
        assert_eq!(result_mapped, expected_mapped);
        // if result_mapped != expected_mapped {
        //     println!(
        //         "original: {}:{}-{} expected: {:?}  actual: {:?}",
        //         expected.original_chrom,
        //         expected.original_start,
        //         expected.original_end,
        //         expected_mapped,
        //         result_mapped
        //     )
        // }
    }

    // TODO: test changes field
    Ok(())
}

#[allow(clippy::unreadable_literal)]
#[test]
fn test_lift_region_changes() -> Result<(), LiftOverError> {
    let lift_over = PositionLiftOver::load(
        &include_bytes!("../../testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain")[..],
    )?;

    let lifted = lift_over.lift_region("chr22", 49149136..49149176);
    assert_eq!(lifted.len(), 1);
    assert_eq!(lifted[0].chromosome.name, "chr22");
    assert_eq!(lifted[0].start, 49544873);
    assert_eq!(lifted[0].end, 49544951);
    assert_eq!(
        lifted[0].changes,
        vec![
            RegionChangeOp::Aligned(20),
            RegionChangeOp::Insertion(38),
            RegionChangeOp::Aligned(20),
        ]
    );

    let lifted = lift_over.lift_region("chr22", 49149659..49149698);
    assert_eq!(lifted.len(), 1);
    assert_eq!(lifted[0].chromosome.name, "chr22");
    assert_eq!(lifted[0].start, 49545434);
    assert_eq!(lifted[0].end, 49545468);
    assert_eq!(
        lifted[0].changes,
        vec![
            RegionChangeOp::Aligned(15),
            RegionChangeOp::Deletion(5),
            RegionChangeOp::Aligned(19),
        ]
    );

    let lifted = lift_over.lift_region("chr22", 49149047..49149781);
    assert_eq!(lifted[0].chromosome.name, "chr22");
    assert_eq!(lifted[0].start, 49544784);
    assert_eq!(lifted[0].end, 49545551);
    assert_eq!(
        lifted[0].changes,
        vec![
            RegionChangeOp::Aligned(109),
            RegionChangeOp::Insertion(38),
            RegionChangeOp::Aligned(518),
            RegionChangeOp::Deletion(5),
            RegionChangeOp::Aligned(102),
        ]
    );
    Ok(())
}
