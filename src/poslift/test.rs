use super::*;

use crate::LiftOverError;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead};

#[test]
fn test_target_region() {
    let value = TargetRegion {
        reference_chromosome_index: 1,
        query_chromosome_index: 1,
        chain_index: 1,
        reference_start: 10,
        reference_end: 20,
        query_start: 30,
        query_end: 40,
        strand: Strand::Forward,
        is_in_gap: true,
    };
    assert_eq!(value.contains_reference_position(9), false);
    assert_eq!(value.contains_reference_position(10), true);
    assert_eq!(value.contains_reference_position(19), true);
    assert_eq!(value.contains_reference_position(20), false);

    assert_eq!(value.contains_query_position(29), false);
    assert_eq!(value.contains_query_position(30), true);
    assert_eq!(value.contains_query_position(39), true);
    assert_eq!(value.contains_query_position(40), false);

    assert_eq!(value.is_indel(), false);
}

#[test]
fn test_target_region2() {
    let value = TargetRegion {
        reference_chromosome_index: 1,
        query_chromosome_index: 1,
        chain_index: 1,
        reference_start: 10,
        reference_end: 20,
        query_start: 30,
        query_end: 50,
        strand: Strand::Forward,
        is_in_gap: true,
    };
    assert_eq!(value.contains_reference_position(9), false);
    assert_eq!(value.contains_reference_position(10), true);
    assert_eq!(value.contains_reference_position(19), true);
    assert_eq!(value.contains_reference_position(20), false);

    assert_eq!(value.contains_query_position(29), false);
    assert_eq!(value.contains_query_position(30), true);
    assert_eq!(value.contains_query_position(40), true);
    assert_eq!(value.contains_query_position(39), true);
    assert_eq!(value.contains_query_position(50), false);

    assert_eq!(value.is_indel(), true);
}

#[allow(clippy::unreadable_literal)]
#[allow(clippy::cognitive_complexity)]
#[test]
fn test_liftover1() -> Result<(), LiftOverError> {
    let hg19_hg38_reader =
        MultiGzDecoder::new(File::open("testfiles/hg19ToHg38/hg19ToHg38.over.chain.gz")?);
    let hg19_hg38_liftover = PositionLiftOver::load(hg19_hg38_reader)?;

    assert_eq!(
        hg19_hg38_liftover.lift_position("chr1", 7424071),
        vec![LiftOverResult {
            chromosome: &hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr1")
                .unwrap(),
            position: 7364011,
            chain_index: 0,
            strand: Strand::Forward,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_position("chrY", 58832218),
        vec![LiftOverResult {
            chromosome: &hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chrY")
                .unwrap(),
            position: 56758651,
            chain_index: 1255,
            strand: Strand::Reverse,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_position("chrY", 58906573),
        vec![LiftOverResult {
            chromosome: &hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chrY")
                .unwrap(),
            position: 56684296,
            chain_index: 1255,
            strand: Strand::Reverse,
        }]
    );

    {
        let search_region_result = hg19_hg38_liftover.search_target("chr8", 1243018, 1243021);
        assert_eq!(search_region_result.len(), 1);
        assert_eq!(search_region_result[0].len(), 3);
        for i in 0..3 {
            assert_eq!(
                hg19_hg38_liftover.reference_chromosomes()
                    [search_region_result[0][i].reference_chromosome_index]
                    .name,
                "chr8"
            );
        }
        assert_eq!(search_region_result[0][1].reference_start, 1243019);
        assert_eq!(search_region_result[0][1].reference_end, 1243020);
        assert_eq!(search_region_result[0][2].reference_start, 1243020);
    }

    assert_eq!(
        hg19_hg38_liftover.lift_position("chr20", 52009099),
        vec![LiftOverResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr20")
                .unwrap(),
            position: 53392561,
            chain_index: 689,
            strand: Strand::Forward,
        }]
    );

    {
        let search_region_result = hg19_hg38_liftover.search_target("chr20", 52009099, 52009100);
        //println!("{:?}", search_region_result);
        assert_eq!(search_region_result.len(), 1);
        assert_eq!(search_region_result[0].len(), 1);
        for i in 0..1 {
            assert_eq!(
                hg19_hg38_liftover.reference_chromosomes()
                    [search_region_result[0][i].reference_chromosome_index]
                    .name,
                "chr20"
            );
        }
        assert_eq!(search_region_result[0][0].reference_start, 52009099);
        assert_eq!(search_region_result[0][0].query_start, 53392561);
        assert_eq!(search_region_result[0][0].is_indel(), false);
    }

    {
        let search_region_result = hg19_hg38_liftover.search_target("chr20", 52009098, 52009099);
        //println!("{:?}", search_region_result);
        assert_eq!(search_region_result.len(), 1);
        assert_eq!(search_region_result[0].len(), 1);
        for i in 0..1 {
            assert_eq!(
                hg19_hg38_liftover.reference_chromosomes()
                    [search_region_result[0][i].reference_chromosome_index]
                    .name,
                "chr20"
            );
        }
        assert_eq!(search_region_result[0][0].reference_end, 52009099);
        assert_eq!(search_region_result[0][0].query_end, 53392560);
        assert_eq!(search_region_result[0][0].is_indel(), false);
    }

    {
        let search_region_result = hg19_hg38_liftover.search_target("chr20", 52009098, 52009100);
        //println!("{:?}", search_region_result);
        assert_eq!(search_region_result.len(), 1);
        assert_eq!(search_region_result[0].len(), 3);
        for i in 0..2 {
            assert_eq!(
                hg19_hg38_liftover.reference_chromosomes()
                    [search_region_result[0][i].reference_chromosome_index]
                    .name,
                "chr20"
            );
        }
        assert_eq!(search_region_result[0][0].reference_end, 52009099);
        assert_eq!(search_region_result[0][1].reference_start, 52009099);
        assert_eq!(search_region_result[0][1].reference_end, 52009099);
        assert_eq!(search_region_result[0][2].reference_start, 52009099);

        assert_eq!(search_region_result[0][0].query_end, 53392560);
        assert_eq!(search_region_result[0][1].query_start, 53392560);
        assert_eq!(search_region_result[0][1].query_end, 53392561);
        assert_eq!(search_region_result[0][2].query_start, 53392561);
    }

    // chr4
    assert_eq!(
        hg19_hg38_liftover.lift_position("chr4_ctg9_hap1", 320130),
        Vec::new()
    );
    assert_eq!(
        hg19_hg38_liftover.lift_position("chr4_ctg9_hap1", 320129),
        vec![LiftOverResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr4")
                .unwrap(),
            position: 68742137,
            chain_index: 758,
            strand: Strand::Forward,
        }]
    );
    assert_eq!(
        hg19_hg38_liftover.lift_position("chr4_ctg9_hap1", 320131),
        vec![LiftOverResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr4")
                .unwrap(),
            position: 68742138,
            chain_index: 758,
            strand: Strand::Forward,
        }]
    );

    {
        let search_region_result =
            hg19_hg38_liftover.search_target("chr4_ctg9_hap1", 320130, 320132);
        //println!("{:?}", search_region_result);
        assert_eq!(search_region_result.len(), 1);
        assert_eq!(search_region_result[0].len(), 2);
        for i in 0..2 {
            assert_eq!(
                hg19_hg38_liftover.reference_chromosomes()
                    [search_region_result[0][i].reference_chromosome_index]
                    .name,
                "chr4_ctg9_hap1"
            );
            assert_eq!(
                hg19_hg38_liftover.query_chromosomes()
                    [search_region_result[0][i].query_chromosome_index]
                    .name,
                "chr4"
            );
        }
        assert_eq!(search_region_result[0][0].reference_start, 320130);
        assert_eq!(search_region_result[0][0].reference_end, 320131);
        assert_eq!(search_region_result[0][1].reference_start, 320131);

        assert_eq!(search_region_result[0][0].query_start, 68742138);
        assert_eq!(search_region_result[0][0].query_end, 68742138);
        assert_eq!(search_region_result[0][1].query_start, 68742138);
    }

    Ok(())
}

#[test]
fn test_liftover2() -> Result<(), LiftOverError> {
    let hg19_hg38_reader =
        MultiGzDecoder::new(File::open("testfiles/hg19ToHg38/hg19ToHg38.over.chain.gz")?);
    let hg19_hg38_liftover = PositionLiftOver::load(hg19_hg38_reader)?;

    let mut test_loader = io::BufReader::new(File::open(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.test.txt",
    )?);
    let mut line = String::new();
    test_loader.read_line(&mut line)?; // skip header

    loop {
        line.clear();
        let read_bytes = test_loader.read_line(&mut line)?;
        if read_bytes == 0 {
            break;
        }
        let elements: Vec<_> = line.split('\t').collect();

        let mut liftover_results: Vec<_> = hg19_hg38_liftover
            .lift_position(elements[0], elements[1].parse::<u64>().unwrap())
            .iter()
            .map(|x| (x.chromosome.name.to_string(), x.position, x.strand))
            .collect();
        liftover_results.sort();

        let mut expected_results: Vec<_> = if elements[2].is_empty() {
            vec![]
        } else {
            elements[2]
                .trim()
                .split(',')
                .zip(
                    elements[3]
                        .trim()
                        .split(',')
                        .zip(elements[4].trim().split(',')),
                )
                .map(|x| {
                    (
                        x.0.to_string(),
                        (x.1).0.parse::<u64>().unwrap(),
                        match (x.1).1 {
                            "+" => Strand::Forward,
                            "-" => Strand::Reverse,
                            _ => {
                                eprintln!("Invalid char: \"{}\"", (x.1).1);
                                unreachable!()
                            }
                        },
                    )
                })
                .collect()
        };
        expected_results.sort();

        assert_eq!(liftover_results, expected_results);
    }

    Ok(())
}

#[allow(clippy::unreadable_literal)]
#[test]
fn test_lift_region() -> Result<(), LiftOverError> {
    let hg19_hg38_reader =
        MultiGzDecoder::new(File::open("testfiles/hg19ToHg38/hg19ToHg38.over.chain.gz")?);
    let hg19_hg38_liftover = PositionLiftOver::load(hg19_hg38_reader)?;

    //single target region
    assert_eq!(
        hg19_hg38_liftover.lift_region("chr20", 52009088, 52009098),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr20")
                .unwrap(),
            start: 53392549,
            end: 53392559,
            strand: Strand::Forward,
            changes: vec![RegionChangeOp::Aligned(10)],
            chain_index: 689,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_region("chr20", 52009088, 52009098),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr20")
                .unwrap(),
            start: 53392549,
            end: 53392559,
            strand: Strand::Forward,
            changes: vec![RegionChangeOp::Aligned(10)],
            chain_index: 689,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_region("chr20", 52009099, 52009109),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr20")
                .unwrap(),
            start: 53392561,
            end: 53392571,
            strand: Strand::Forward,
            changes: vec![RegionChangeOp::Aligned(10)],
            chain_index: 689,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_region("chr20", 52009100, 52009110),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr20")
                .unwrap(),
            start: 53392562,
            end: 53392572,
            strand: Strand::Forward,
            changes: vec![RegionChangeOp::Aligned(10)],
            chain_index: 689,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_region("chr22", 16071914, 16071924),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr22")
                .unwrap(),
            start: 15906030,
            end: 15906040,
            strand: Strand::Reverse,
            changes: vec![RegionChangeOp::Aligned(10)],
            chain_index: 696,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_region("chr1", 12812500, 12812506),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr1")
                .unwrap(),
            start: 12752552,
            end: 12752557,
            strand: Strand::Forward,
            changes: vec![
                RegionChangeOp::Aligned(3),
                RegionChangeOp::Deletion(1),
                RegionChangeOp::Aligned(2)
            ],
            chain_index: 0,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_region("chr1", 12803239, 12803244),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr1")
                .unwrap(),
            start: 12743290,
            end: 12743296,
            strand: Strand::Forward,
            changes: vec![
                RegionChangeOp::Aligned(2),
                RegionChangeOp::Insertion(1),
                RegionChangeOp::Aligned(3)
            ],
            chain_index: 0,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_region("chr10", 116895603, 116904040),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr10")
                .unwrap(),
            start: 115135840,
            end: 115144276,
            strand: Strand::Forward,
            changes: vec![
                RegionChangeOp::Aligned(95),
                RegionChangeOp::Deletion(1),
                RegionChangeOp::Aligned(19),
                RegionChangeOp::Insertion(1),
                RegionChangeOp::Deletion(1),
                RegionChangeOp::Aligned(7456),
                RegionChangeOp::Insertion(1),
                RegionChangeOp::Aligned(571),
                RegionChangeOp::Deletion(1),
                RegionChangeOp::Aligned(293),
            ],
            chain_index: 184,
        }]
    );

    assert_eq!(
        hg19_hg38_liftover.lift_region("chr22", 16067539, 16072072),
        vec![LiftRegionResult {
            chromosome: hg19_hg38_liftover
                .chain_file
                .query_chromosome_by_name("chr22")
                .unwrap(),
            start: 15905890,
            end: 15910417,
            strand: Strand::Reverse,
            changes: vec![
                RegionChangeOp::Aligned(153),
                RegionChangeOp::Insertion(1),
                RegionChangeOp::Deletion(1),
                RegionChangeOp::Aligned(25),
                RegionChangeOp::Insertion(1),
                RegionChangeOp::Deletion(1),
                RegionChangeOp::Aligned(2051),
                RegionChangeOp::Insertion(14),
                RegionChangeOp::Deletion(13),
                RegionChangeOp::Aligned(211),
                RegionChangeOp::Insertion(1),
                RegionChangeOp::Aligned(1931),
                RegionChangeOp::Deletion(8),
                RegionChangeOp::Aligned(139),
            ],
            chain_index: 696,
        }]
    );

    Ok(())
}

#[test]
fn test_lift2() -> Result<(), LiftOverError> {
    let lift = PositionLiftOver::load(File::open(
        "testfiles/lift-variant-test/lifttest/seq-a--to--seq-b.chain",
    )?)?;

    let mut lifted = lift.lift_region("seq-a", 11540, 11543);
    lifted.sort_by_key(|x| x.start);
    assert_eq!(
        lifted,
        vec![
            LiftRegionResult {
                chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
                start: 7400,
                end: 7403,
                strand: Strand::Forward,
                changes: vec![RegionChangeOp::Aligned(3)],
                chain_index: 2,
            },
            LiftRegionResult {
                chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
                start: 15915,
                end: 15918,
                strand: Strand::Forward,
                changes: vec![RegionChangeOp::Aligned(3)],
                chain_index: 3,
            },
        ]
    );

    let mut lifted = lift.lift_region("seq-a", 11542, 11549);
    lifted.sort_by_key(|x| x.start);
    assert_eq!(
        lifted,
        vec![
            LiftRegionResult {
                chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
                start: 7402,
                end: 7407,
                strand: Strand::Forward,
                changes: vec![
                    RegionChangeOp::Aligned(2),
                    RegionChangeOp::Deletion(2),
                    RegionChangeOp::Aligned(3)
                ],
                chain_index: 2,
            },
            LiftRegionResult {
                chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
                start: 15917,
                end: 15924,
                strand: Strand::Forward,
                changes: vec![RegionChangeOp::Aligned(7)],
                chain_index: 3,
            },
        ]
    );

    let mut lifted = lift.lift_region("seq-a", 11542, 11547);
    lifted.sort_by_key(|x| x.start);
    assert_eq!(
        lifted[0],
        LiftRegionResult {
            chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
            start: 7402,
            end: 7405,
            strand: Strand::Forward,
            changes: vec![
                RegionChangeOp::Aligned(2),
                RegionChangeOp::Deletion(2),
                RegionChangeOp::Aligned(1)
            ],
            chain_index: 2,
        },
    );

    let mut lifted = lift.lift_region("seq-a", 11542, 11546);
    lifted.sort_by_key(|x| x.start);
    assert_eq!(
        lifted[0],
        LiftRegionResult {
            chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
            start: 7402,
            end: 7404,
            strand: Strand::Forward,
            changes: vec![RegionChangeOp::Aligned(2), RegionChangeOp::Deletion(2),],
            chain_index: 2,
        },
    );

    let mut lifted = lift.lift_region("seq-a", 11542, 11545);
    lifted.sort_by_key(|x| x.start);
    assert_eq!(
        lifted[0],
        LiftRegionResult {
            chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
            start: 7402,
            end: 7404,
            strand: Strand::Forward,
            changes: vec![RegionChangeOp::Aligned(2), RegionChangeOp::Deletion(1),],
            chain_index: 2,
        },
    );

    let mut lifted = lift.lift_region("seq-a", 11542, 11544);
    lifted.sort_by_key(|x| x.start);
    assert_eq!(
        lifted[0],
        LiftRegionResult {
            chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
            start: 7402,
            end: 7404,
            strand: Strand::Forward,
            changes: vec![RegionChangeOp::Aligned(2),],
            chain_index: 2,
        },
    );

    let mut lifted = lift.lift_region("seq-a", 11572, 11576);
    lifted.sort_by_key(|x| x.start);
    assert_eq!(
        lifted,
        vec![
            LiftRegionResult {
                chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
                start: 7430,
                end: 7436,
                strand: Strand::Forward,
                changes: vec![
                    RegionChangeOp::Aligned(2),
                    RegionChangeOp::Insertion(2),
                    RegionChangeOp::Aligned(2),
                ],
                chain_index: 2,
            },
            LiftRegionResult {
                chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
                start: 15947,
                end: 15951,
                strand: Strand::Forward,
                changes: vec![RegionChangeOp::Aligned(4),],
                chain_index: 3,
            },
        ]
    );

    assert_eq!(
        lift.lift_region("seq-a", 31669, 31673),
        vec![LiftRegionResult {
            chromosome: lift.chain_file.query_chromosome_by_name("seq-b").unwrap(),
            start: 33514,
            end: 33518,
            strand: Strand::Reverse,
            changes: vec![RegionChangeOp::Aligned(4),],
            chain_index: 0,
        },]
    );

    // TODO: Add more edge cases

    //TODO: check here
    // assert_eq!(
    //     lift.lift_region("seq-a", 31692, 31700),
    //     vec![LiftRegionResult {
    //         chromosome: lift.chain_file.query_chromosome_by_name("seq-b"),
    //         start: 33487,
    //         end: 33492,
    //         strand: Strand::Reverse,
    //         changes: vec![
    //             RegionChangeOp::Aligned(2),
    //             RegionChangeOp::Deletion(3),
    //             RegionChangeOp::Aligned(3)
    //         ],
    //         chain_index: 0,
    //     },]
    // );

    Ok(())
}
