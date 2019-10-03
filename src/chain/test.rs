use super::*;

use bio::io::fasta::IndexedReader;
use flate2::read::MultiGzDecoder;
use std::fs::File;

#[allow(clippy::unreadable_literal)]
#[test]
fn test_load() -> Result<(), LiftOverError> {
    let hg19_hg38_reader =
        MultiGzDecoder::new(File::open("testfiles/hg19ToHg38/hg19ToHg38.over.chain.gz")?);
    let hg19_hg38_liftover = ChainFile::load(hg19_hg38_reader)?;

    assert_eq!(hg19_hg38_liftover.chain_list.len(), 1278);

    assert_eq!(hg19_hg38_liftover.chain_list[0].chain_interval.len(), 4115);
    assert_eq!(hg19_hg38_liftover.chain_list[1].chain_interval.len(), 1459);
    assert_eq!(hg19_hg38_liftover.chain_list[1277].chain_interval.len(), 2);

    // check first chain
    let chain_without_interval_list = Chain {
        chain_interval: Vec::new(),
        ..hg19_hg38_liftover.chain_list[0].clone()
    };

    assert_eq!(
        chain_without_interval_list,
        Chain {
            chain_interval: Vec::new(),
            score: 20851231461,
            reference_chromosome: Chromosome {
                name: "chr1".to_string(),
                length: 249250621,
            },
            reference_strand: Strand::Forward,
            reference_start: 10000,
            reference_end: 249240621,
            query_chromosome: Chromosome {
                name: "chr1".to_string(),
                length: 248956422,
            },
            query_strand: Strand::Forward,
            query_start: 10000,
            query_end: 248946422,
            chain_id: "2".to_string(),
        }
    );

    assert_eq!(
        hg19_hg38_liftover.chain_list[0].chain_interval[0],
        ChainInterval {
            size: 167376,
            difference_reference: Some(50041),
            difference_query: Some(80290),
        }
    );

    // check second chain
    let chain_without_interval_list = Chain {
        chain_interval: Vec::new(),
        ..hg19_hg38_liftover.chain_list[1].clone()
    };

    assert_eq!(
        chain_without_interval_list,
        Chain {
            chain_interval: Vec::new(),
            score: 51897113,
            reference_chromosome: Chromosome {
                name: "chr1".to_string(),
                length: 249250621,
            },
            reference_strand: Strand::Forward,
            reference_start: 144274511,
            reference_end: 149034137,
            query_chromosome: Chromosome {
                name: "chr1".to_string(),
                length: 248956422,
            },
            query_strand: Strand::Reverse,
            query_start: 99285234,
            query_end: 105680422,
            chain_id: "98".to_string(),
        }
    );

    // check 6th chain
    assert_eq!(
        hg19_hg38_liftover.chain_list[5].clone(),
        Chain {
            chain_interval: vec![ChainInterval {
                size: 186739,
                difference_reference: None,
                difference_query: None,
            }],
            score: 17644093,
            reference_chromosome: Chromosome {
                name: "chr1".to_string(),
                length: 249250621,
            },
            reference_strand: Strand::Forward,
            reference_start: 142781022,
            reference_end: 142967761,
            query_chromosome: Chromosome {
                name: "chrUn_KI270742v1".to_string(),
                length: 186739,
            },
            query_strand: Strand::Forward,
            query_start: 0,
            query_end: 186739,
            chain_id: "358".to_string(),
        }
    );

    // check last chain
    assert_eq!(
        hg19_hg38_liftover.chain_list[1277].clone(),
        Chain {
            chain_interval: vec![
                ChainInterval {
                    size: 41,
                    difference_reference: Some(1338),
                    difference_query: Some(158),
                },
                ChainInterval {
                    size: 27,
                    difference_reference: None,
                    difference_query: None,
                }
            ],
            score: 945,
            reference_chromosome: Chromosome {
                name: "chrY".to_string(),
                length: 59373566,
            },
            reference_strand: Strand::Forward,
            reference_start: 1135403,
            reference_end: 1136809,
            query_chromosome: Chromosome {
                name: "chrX".to_string(),
                length: 156040895,
            },
            query_strand: Strand::Forward,
            query_start: 1085840,
            query_end: 1086066,
            chain_id: "13541154".to_string(),
        }
    );

    Ok(())
}

#[test]
fn test_normalize_chain_file() -> Result<(), LiftOverError> {
    let test_chain = ChainFile::load(File::open(
        "testfiles/lift-variant-test/normalize-test/before-normalized.chain",
    )?)?;
    let mut reference = IndexedReader::from_file(&"testfiles/lift-variant-test/sequence/seq-a.fa")?;
    let mut query = IndexedReader::from_file(&"testfiles/lift-variant-test/sequence/seq-b.fa")?;

    let mut test_writer = File::create("target/write-test.chain")?;
    test_chain.write(&mut test_writer)?;

    let left_aligned_chain = test_chain.left_align(&mut reference, &mut query)?;
    let mut test_writer = File::create("target/write-test2.chain")?;
    left_aligned_chain.write(&mut test_writer)?;

    assert_eq!(
        left_aligned_chain.chain_list[0].chain_interval,
        vec![
            ChainInterval {
                size: 5,
                difference_reference: Some(3),
                difference_query: Some(0),
            },
            ChainInterval {
                size: 13,
                difference_reference: Some(0),
                difference_query: Some(3),
            },
            ChainInterval {
                size: 6,
                difference_reference: None,
                difference_query: None,
            }
        ]
    );

    assert_eq!(
        left_aligned_chain.chain_list[1].chain_interval,
        vec![
            ChainInterval {
                size: 7,
                difference_reference: Some(1),
                difference_query: Some(0)
            },
            ChainInterval {
                size: 11,
                difference_reference: Some(0),
                difference_query: Some(1)
            },
            ChainInterval {
                size: 21,
                difference_reference: None,
                difference_query: None
            },
        ]
    );

    Ok(())
}

#[test]
fn test_chain_interval_valid() {
    assert_eq!(
        ChainInterval {
            size: 0,
            difference_reference: Some(10),
            difference_query: Some(0),
        }
        .is_valid(),
        false
    );

    assert_eq!(
        ChainInterval {
            size: 10,
            difference_reference: Some(0),
            difference_query: Some(0),
        }
        .is_valid(),
        false
    );

    assert_eq!(
        ChainInterval {
            size: 10,
            difference_reference: Some(10),
            difference_query: Some(0),
        }
        .is_valid(),
        true
    );

    assert_eq!(
        ChainInterval {
            size: 10,
            difference_reference: Some(0),
            difference_query: Some(10),
        }
        .is_valid(),
        true
    );

    assert_eq!(
        ChainInterval {
            size: 10,
            difference_reference: Some(10),
            difference_query: Some(10),
        }
        .is_valid(),
        true
    );
}

#[test]
fn test_chain_cleanup_no_change() {
    let chain_intervals = vec![
        ChainInterval {
            size: 10,
            difference_reference: Some(5),
            difference_query: Some(0),
        },
        ChainInterval {
            size: 10,
            difference_reference: Some(5),
            difference_query: Some(10),
        },
        ChainInterval {
            size: 10,
            difference_reference: Some(0),
            difference_query: Some(5),
        },
        ChainInterval {
            size: 10,
            difference_reference: Some(5),
            difference_query: Some(0),
        },
    ];
    let cleanup_chain = Chain::cleanup(&chain_intervals);

    assert_eq!(chain_intervals, cleanup_chain);
}

#[test]
fn test_chain_cleanup_with_change() {
    let chain_intervals = vec![
        ChainInterval {
            size: 4,
            difference_reference: Some(5),
            difference_query: Some(0),
        },
        ChainInterval {
            size: 0,
            difference_reference: Some(5),
            difference_query: Some(13),
        },
        ChainInterval {
            size: 17,
            difference_reference: Some(0),
            difference_query: Some(0),
        },
        ChainInterval {
            size: 19,
            difference_reference: Some(8),
            difference_query: Some(0),
        },
    ];
    let cleanup_chain = Chain::cleanup(&chain_intervals);

    assert_eq!(
        vec![
            ChainInterval {
                size: 4,
                difference_reference: Some(10),
                difference_query: Some(13),
            },
            ChainInterval {
                size: 36,
                difference_reference: Some(8),
                difference_query: Some(0),
            },
        ],
        cleanup_chain
    );
}
