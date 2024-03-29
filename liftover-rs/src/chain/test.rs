use super::*;

use bio::io::fasta::IndexedReader;
use std::fs::File;
//use std::str;

#[allow(clippy::unreadable_literal)]
#[test]
fn test_load() -> Result<(), LiftOverError> {
    let chain_file = ChainFile::load(
        &include_bytes!("../../testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain")[..],
    )?;
    assert_eq!(chain_file.chain_list.len(), 104);

    assert_eq!(chain_file.chain_list[0].chain_interval.len(), 212);
    assert_eq!(chain_file.chain_list[1].chain_interval.len(), 50);
    assert_eq!(chain_file.chain_list[103].chain_interval.len(), 10);

    // check last chain
    assert_eq!(
        chain_file.chain_list[103].clone(),
        Chain {
            chain_interval: vec![
                ChainInterval {
                    size: 1058,
                    difference_original: Some(0),
                    difference_new: Some(1),
                },
                ChainInterval {
                    size: 260,
                    difference_original: Some(0),
                    difference_new: Some(5),
                },
                ChainInterval {
                    size: 384,
                    difference_original: Some(0),
                    difference_new: Some(2),
                },
                ChainInterval {
                    size: 95,
                    difference_original: Some(3),
                    difference_new: Some(0),
                },
                ChainInterval {
                    size: 28,
                    difference_original: Some(0),
                    difference_new: Some(1),
                },
                ChainInterval {
                    size: 1595,
                    difference_original: Some(0),
                    difference_new: Some(7),
                },
                ChainInterval {
                    size: 259,
                    difference_original: Some(5),
                    difference_new: Some(0),
                },
                ChainInterval {
                    size: 140,
                    difference_original: Some(1),
                    difference_new: Some(0),
                },
                ChainInterval {
                    size: 5,
                    difference_original: Some(1),
                    difference_new: Some(0),
                },
                ChainInterval {
                    size: 698,
                    difference_original: None,
                    difference_new: None,
                }
            ],
            score: 4900,
            original_chromosome: Chromosome {
                name: "chr22".to_string(),
                length: 50818468,
            },
            original_strand: Strand::Forward,
            original_start: 12802792,
            original_end: 12807324,
            new_chromosome: Chromosome {
                name: "chr22".to_string(),
                length: 51304566,
            },
            new_strand: Strand::Reverse,
            new_start: 28668699,
            new_end: 28673237,
            chain_id: "104".to_string(),
        }
    );

    Ok(())
}

#[test]
fn test_normalize_chain_file() -> anyhow::Result<()> {
    let test_chain = ChainFile::load(&include_bytes!("before-normalize.chain")[..])?;
    let mut reference =
        IndexedReader::from_file(&"testfiles/genomes/GRCh38/GRCh38.chr22.genome.fa")?;
    let mut query = IndexedReader::from_file(&"testfiles/genomes/GRCh37/GRCh37.chr22.genome.fa")?;

    std::fs::create_dir_all("../target/test-output/")?;

    let mut test_writer = File::create("../target/test-output/write-test.chain")?;
    test_chain.write(&mut test_writer)?;

    let left_aligned_chain = test_chain.left_align(&mut reference, &mut query)?;
    let mut test_writer = File::create("../target/test-output/write-test2.chain")?;
    left_aligned_chain.write(&mut test_writer)?;

    assert_eq!(
        left_aligned_chain.chain_list[0].chain_interval,
        vec![
            ChainInterval {
                size: 3,
                difference_original: Some(2),
                difference_new: Some(0),
            },
            ChainInterval {
                size: 2,
                difference_original: None,
                difference_new: None,
            }
        ]
    );

    assert_eq!(
        left_aligned_chain.chain_list[1].chain_interval,
        vec![
            ChainInterval {
                size: 2,
                difference_original: Some(0),
                difference_new: Some(3),
            },
            ChainInterval {
                size: 5,
                difference_original: None,
                difference_new: None,
            }
        ]
    );

    assert_eq!(
        left_aligned_chain.chain_list[2].chain_interval,
        vec![
            ChainInterval {
                size: 2,
                difference_original: Some(0),
                difference_new: Some(3),
            },
            ChainInterval {
                size: 5,
                difference_original: None,
                difference_new: None,
            }
        ]
    );

    Ok(())
}

#[test]
fn test_chain_interval_valid() {
    assert_eq!(
        ChainInterval {
            size: 0,
            difference_original: Some(10),
            difference_new: Some(0),
        }
        .is_valid(),
        false
    );

    assert_eq!(
        ChainInterval {
            size: 10,
            difference_original: Some(0),
            difference_new: Some(0),
        }
        .is_valid(),
        false
    );

    assert_eq!(
        ChainInterval {
            size: 10,
            difference_original: Some(10),
            difference_new: Some(0),
        }
        .is_valid(),
        true
    );

    assert_eq!(
        ChainInterval {
            size: 10,
            difference_original: Some(0),
            difference_new: Some(10),
        }
        .is_valid(),
        true
    );

    assert_eq!(
        ChainInterval {
            size: 10,
            difference_original: Some(10),
            difference_new: Some(10),
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
            difference_original: Some(5),
            difference_new: Some(0),
        },
        ChainInterval {
            size: 10,
            difference_original: Some(5),
            difference_new: Some(10),
        },
        ChainInterval {
            size: 10,
            difference_original: Some(0),
            difference_new: Some(5),
        },
        ChainInterval {
            size: 10,
            difference_original: Some(5),
            difference_new: Some(0),
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
            difference_original: Some(5),
            difference_new: Some(0),
        },
        ChainInterval {
            size: 0,
            difference_original: Some(5),
            difference_new: Some(13),
        },
        ChainInterval {
            size: 17,
            difference_original: Some(0),
            difference_new: Some(0),
        },
        ChainInterval {
            size: 19,
            difference_original: Some(8),
            difference_new: Some(0),
        },
    ];
    let cleanup_chain = Chain::cleanup(&chain_intervals);

    assert_eq!(
        vec![
            ChainInterval {
                size: 4,
                difference_original: Some(10),
                difference_new: Some(13),
            },
            ChainInterval {
                size: 36,
                difference_original: Some(8),
                difference_new: Some(0),
            },
        ],
        cleanup_chain
    );
}
