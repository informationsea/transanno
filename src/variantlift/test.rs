use super::*;
use std::fs::File;

#[test]
fn test_lift_variant() -> Result<(), LiftOverError> {
    let reference_fasta =
        bio::io::fasta::IndexedReader::from_file(&"testfiles/lift-variant-test/sequence/seq-a.fa")?;
    let query_fasta =
        bio::io::fasta::IndexedReader::from_file(&"testfiles/lift-variant-test/sequence/seq-b.fa")?;

    let mut variant_liftover = VariantLiftOver::load(
        File::open("testfiles/lift-variant-test/lifttest/seq-a--to--seq-b.chain")?,
        reference_fasta,
        query_fasta,
    )?;

    assert_eq!(
        variant_liftover
            .lift_variant(
                &Variant {
                    chromosome: "seq-a".to_string(),
                    position: 9904,
                    reference: b"C".to_vec(),
                    alternative: vec![b"A".to_vec()]
                },
                100,
                100
            )
            .unwrap(),
        vec![Err(
            error::VariantLiftOverError::ReferenceSequenceIsNotMatch
        )]
    );

    let mut lifted = variant_liftover.lift_variant(
        &Variant {
            chromosome: "seq-a".to_string(),
            position: 11542,
            reference: b"C".to_vec(),
            alternative: vec![b"T".to_vec()],
        },
        100,
        100,
    )?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 7402,
                strand: Strand::Forward,
                original_reference: b"C".to_vec(),
                reference: b"C".to_vec(),
                alternative: vec![b"T".to_vec()],
                reference_changed: false,
            }),
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 15917,
                strand: Strand::Forward,
                original_reference: b"C".to_vec(),
                reference: b"C".to_vec(),
                alternative: vec![b"T".to_vec()],
                reference_changed: false,
            }),
        ]
    );

    let mut lifted = variant_liftover.lift_variant(
        &Variant {
            chromosome: "seq-a".to_string(),
            position: 11543,
            reference: b"TA".to_vec(),
            alternative: vec![b"T".to_vec(), b"*".to_vec()],
        },
        100,
        100,
    )?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 7403,
                strand: Strand::Forward,
                original_reference: b"TAA".to_vec(),
                reference: b"T".to_vec(),
                alternative: vec![b"TA".to_vec(), b"*".to_vec()],
                reference_changed: true,
            }),
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 15918,
                strand: Strand::Forward,
                original_reference: b"TA".to_vec(),
                reference: b"TA".to_vec(),
                alternative: vec![b"T".to_vec(), b"*".to_vec()],
                reference_changed: false,
            }),
        ]
    );

    let mut lifted =
        variant_liftover.lift_variant(&Variant::new("seq-a", 11543, b"T", &[b"TA"]), 100, 100)?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 7403,
                strand: Strand::Forward,
                original_reference: b"TAA".to_vec(),
                reference: b"T".to_vec(),
                alternative: vec![b"TAAA".to_vec()],
                reference_changed: true,
            }),
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 15918,
                strand: Strand::Forward,
                original_reference: b"T".to_vec(),
                reference: b"T".to_vec(),
                alternative: vec![b"TA".to_vec()],
                reference_changed: false,
            }),
        ]
    );

    let mut lifted =
        variant_liftover.lift_variant(&Variant::new("seq-a", 11573, b"TA", &[b"T"]), 100, 100)?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 7431,
                strand: Strand::Forward,
                original_reference: b"TA".to_vec(),
                reference: b"TAAA".to_vec(),
                alternative: vec![b"T".to_vec()],
                reference_changed: true,
            }),
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 15948,
                strand: Strand::Forward,
                original_reference: b"TA".to_vec(),
                reference: b"TA".to_vec(),
                alternative: vec![b"T".to_vec()],
                reference_changed: false,
            }),
        ]
    );

    let mut lifted =
        variant_liftover.lift_variant(&Variant::new("seq-a", 11573, b"T", &[b"TA"]), 100, 100)?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 7431,
                strand: Strand::Forward,
                original_reference: b"T".to_vec(),
                reference: b"TAA".to_vec(),
                alternative: vec![b"TA".to_vec()],
                reference_changed: true,
            }),
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 15948,
                strand: Strand::Forward,
                original_reference: b"T".to_vec(),
                reference: b"T".to_vec(),
                alternative: vec![b"TA".to_vec()],
                reference_changed: false,
            }),
        ]
    );

    let mut lifted =
        variant_liftover.lift_variant(&Variant::new("seq-a", 15516, b"T", &[b"G"]), 100, 100)?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 11376,
                strand: Strand::Forward,
                original_reference: b"T".to_vec(),
                reference: b"T".to_vec(),
                alternative: vec![b"G".to_vec()],
                reference_changed: false,
            }),
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 19891,
                strand: Strand::Forward,
                original_reference: b"T".to_vec(),
                reference: b"G".to_vec(),
                alternative: vec![b"G".to_vec()],
                reference_changed: true,
            }),
        ]
    );

    let mut lifted = variant_liftover.lift_variant(
        &Variant::new("seq-a", 16396, b"G", &[b"T", b"GA"]),
        100,
        100,
    )?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 12256,
                strand: Strand::Forward,
                original_reference: b"G".to_vec(),
                reference: b"G".to_vec(),
                alternative: vec![b"T".to_vec(), b"GA".to_vec()],
                reference_changed: false,
            }),
            Ok(LiftedVariant {
                chromosome: "seq-b".to_string(),
                position: 20769,
                strand: Strand::Forward,
                original_reference: b"CAGA".to_vec(),
                reference: b"C".to_vec(),
                alternative: vec![b"CATA".to_vec(), b"CAGAA".to_vec()],
                reference_changed: true,
            }),
        ]
    );

    Ok(())
}

#[test]
fn test_lift_variant_reverse() -> Result<(), LiftOverError> {
    let mut reference_fasta =
        bio::io::fasta::IndexedReader::from_file(&"testfiles/lift-variant-test/sequence/seq-a.fa")?;
    let mut query_fasta =
        bio::io::fasta::IndexedReader::from_file(&"testfiles/lift-variant-test/sequence/seq-b.fa")?;
    let chain_file = ChainFile::load(File::open(
        "testfiles/lift-variant-test/lifttest/seq-a--to--seq-b.chain",
    )?)?;

    let mut variant_liftover = VariantLiftOver::new(
        chain_file.left_align(&mut reference_fasta, &mut query_fasta)?,
        reference_fasta,
        query_fasta,
    );

    let mut lifted =
        variant_liftover.lift_variant(&Variant::new("seq-a", 24384, b"A", &[b"G"]), 100, 100)?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![Ok(LiftedVariant {
            chromosome: "seq-b".to_string(),
            position: 28816,
            strand: Strand::Reverse,
            original_reference: b"T".to_vec(),
            reference: b"T".to_vec(),
            alternative: vec![b"C".to_vec()],
            reference_changed: false,
        }),]
    );

    let mut lifted =
        variant_liftover.lift_variant(&Variant::new("seq-a", 24386, b"T", &[b"G"]), 100, 100)?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![Ok(LiftedVariant {
            chromosome: "seq-b".to_string(),
            position: 28814,
            strand: Strand::Reverse,
            original_reference: b"A".to_vec(),
            reference: b"C".to_vec(),
            alternative: vec![b"C".to_vec()],
            reference_changed: true,
        }),]
    );

    // Add INDEL
    let mut lifted = variant_liftover.lift_variant(
        &Variant::new("seq-a", 25010, b"T", &[b"TA", b"TAA", b"TAAA", b"*"]),
        100,
        100,
    )?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![Ok(LiftedVariant {
            chromosome: "seq-b".to_string(),
            position: 28189,
            strand: Strand::Reverse,
            original_reference: b"G".to_vec(),
            reference: b"G".to_vec(),
            alternative: vec![
                b"GT".to_vec(),
                b"GTT".to_vec(),
                b"GTTT".to_vec(),
                b"*".to_vec()
            ],
            reference_changed: false,
        }),]
    );

    let mut lifted = variant_liftover.lift_variant(
        &Variant::new("seq-a", 25019, b"G", &[b"GA", b"GAA", b"GAAA"]),
        100,
        100,
    )?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![Ok(LiftedVariant {
            chromosome: "seq-b".to_string(),
            position: 28174,
            strand: Strand::Reverse,
            original_reference: b"G".to_vec(),
            reference: b"GTTT".to_vec(),
            alternative: vec![b"GT".to_vec(), b"GTT".to_vec(), b"GTTT".to_vec()],
            reference_changed: true,
        }),]
    );

    let mut lifted = variant_liftover.lift_variant(
        &Variant::new("seq-a", 25038, b"T", &[b"TG", b"C"]),
        100,
        100,
    )?;
    lifted.sort();
    assert_eq!(
        lifted,
        vec![Ok(LiftedVariant {
            chromosome: "seq-b".to_string(),
            position: 28159,
            strand: Strand::Reverse,
            original_reference: b"AAAA".to_vec(),
            reference: b"A".to_vec(),
            alternative: vec![b"AAACA".to_vec(), b"AAAG".to_vec()],
            reference_changed: true,
        }),]
    );

    Ok(())
}

#[allow(clippy::unreadable_literal)]
#[test]
fn test_in_real_data() -> Result<(), LiftOverError> {
    let mut reference_fasta =
        bio::io::fasta::IndexedReader::from_file(&"testfiles/genome/hg19-chr22.fa")?;
    let mut query_fasta =
        bio::io::fasta::IndexedReader::from_file(&"testfiles/genome/hg38-chr22.fa")?;
    let chain_file = ChainFile::load(File::open(
        "testfiles/hg19ToHg38/hg19ToHg38.over.chain.chr22",
    )?)?;

    let mut variant_liftover = VariantLiftOver::new(
        chain_file.left_align(&mut reference_fasta, &mut query_fasta)?,
        reference_fasta,
        query_fasta,
    );

    assert_eq!(
        variant_liftover.lift_variant(&Variant::new("chr22", 17119486, b"G", &[b"C"]), 100, 100)?,
        vec![Ok(LiftedVariant {
            chromosome: "chr22".to_string(),
            position: 16638596,
            strand: Strand::Forward,
            original_reference: b"G".to_vec(),
            reference: b"G".to_vec(),
            alternative: vec![b"C".to_vec()],
            reference_changed: false,
        }),]
    );

    assert_eq!(
        variant_liftover.lift_variant(&Variant::new("chr22", 34955688, b"G", &[b"T"]), 100, 100)?,
        vec![Ok(LiftedVariant {
            chromosome: "chr22".to_string(),
            position: 34559696,
            strand: Strand::Forward,
            original_reference: b"G".to_vec(),
            reference: b"G".to_vec(),
            alternative: vec![b"T".to_vec()],
            reference_changed: false,
        }),]
    );

    assert_eq!(
        variant_liftover.lift_variant(&Variant::new("chr22", 16473784, b"A", &[b"C"]), 100, 100)?,
        vec![Ok(LiftedVariant {
            chromosome: "chr22".to_string(),
            position: 15504177,
            strand: Strand::Reverse,
            original_reference: b"T".to_vec(),
            reference: b"T".to_vec(),
            alternative: vec![b"G".to_vec()],
            reference_changed: false,
        }),]
    );

    assert_eq!(
        variant_liftover.lift_variant(&Variant::new("chr22", 16581406, b"T", &[b"C"]), 100, 100)?,
        vec![Ok(LiftedVariant {
            chromosome: "chr22".to_string(),
            position: 15396555,
            strand: Strand::Reverse,
            original_reference: b"A".to_vec(),
            reference: b"A".to_vec(),
            alternative: vec![b"G".to_vec()],
            reference_changed: false,
        }),]
    );

    // assert_eq!(
    //     variant_liftover.lift_variant(&Variant::new("chr22", 42522750, b"C", &[b"T"]), 3, 3)?,
    //     vec![Ok(LiftedVariant {
    //         chromosome: "chr22".to_string(),
    //         position: 42126748,
    //         strand: Strand::Forward,
    //         original_reference: b"C".to_vec(),
    //         reference: b"C".to_vec(),
    //         alternative: vec![b"T".to_vec()],
    //         reference_changed: false,
    //     }),]
    // );

    Ok(())
}
