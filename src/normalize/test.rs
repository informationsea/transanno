use super::*;
use crate::LiftOverError;
use bio::io::fasta::IndexedReader;

#[test]
fn test_extend_to_left_if_empty_allele_exists() -> Result<(), LiftOverError> {
    let mut fasta = IndexedReader::from_file(&"testfiles/lift-variant-test/sequence/seq-a.fa")?;
    {
        let chromosome = "seq-a";
        let mut position = 18041;
        let mut alleles = vec![b"".to_vec(), b"T".to_vec(), b"TT".to_vec()];
        assert!(extend_to_left_if_empty_allele_exists(
            chromosome,
            &mut fasta,
            &mut position,
            &mut alleles
        )?);
        assert_eq!(position, 18040);
        assert_eq!(
            alleles,
            vec![b"G".to_vec(), b"GT".to_vec(), b"GTT".to_vec()]
        );
    }

    {
        let chromosome = "seq-a";
        let mut position = 18041;
        let mut alleles = vec![b"".to_vec(), b"T".to_vec(), b"TT".to_vec(), b"*".to_vec()];
        assert!(extend_to_left_if_empty_allele_exists(
            chromosome,
            &mut fasta,
            &mut position,
            &mut alleles
        )?);
        assert_eq!(position, 18040);
        assert_eq!(
            alleles,
            vec![
                b"G".to_vec(),
                b"GT".to_vec(),
                b"GTT".to_vec(),
                b"*".to_vec()
            ]
        );
    }

    {
        let chromosome = "seq-a";
        let mut position = 18040;
        let mut alleles = vec![b"G".to_vec(), b"GT".to_vec(), b"GTT".to_vec()];
        assert!(!extend_to_left_if_empty_allele_exists(
            chromosome,
            &mut fasta,
            &mut position,
            &mut alleles
        )?);
        assert_eq!(position, 18040);
        assert_eq!(
            alleles,
            vec![b"G".to_vec(), b"GT".to_vec(), b"GTT".to_vec()]
        );
    }

    Ok(())
}

#[test]
fn test_truncate_right_most_nucleotide_if_allele_ends_with_same() {
    {
        let mut alleles = vec![b"G".to_vec(), b"GT".to_vec(), b"GTT".to_vec()];
        assert!(!truncate_right_most_nucleotide_if_allele_ends_with_same(
            10,
            &mut alleles
        ));
    }

    {
        let mut alleles = vec![b"GTTT".to_vec(), b"GTTTT".to_vec(), b"GTT".to_vec()];
        assert!(truncate_right_most_nucleotide_if_allele_ends_with_same(
            10,
            &mut alleles
        ));
        assert_eq!(
            alleles,
            vec![b"GT".to_vec(), b"GTT".to_vec(), b"G".to_vec()]
        );
    }

    {
        let mut alleles = vec![
            b"GTTT".to_vec(),
            b"GTTTT".to_vec(),
            b"GTT".to_vec(),
            b"*".to_vec(),
        ];
        assert!(truncate_right_most_nucleotide_if_allele_ends_with_same(
            10,
            &mut alleles
        ));
        assert_eq!(
            alleles,
            vec![
                b"GT".to_vec(),
                b"GTT".to_vec(),
                b"G".to_vec(),
                b"*".to_vec()
            ]
        );
    }

    {
        let mut alleles = vec![b"GTA".to_vec(), b"GTC".to_vec(), b"GTG".to_vec()];
        assert!(!truncate_right_most_nucleotide_if_allele_ends_with_same(
            10,
            &mut alleles
        ));
        assert_eq!(
            alleles,
            vec![b"GTA".to_vec(), b"GTC".to_vec(), b"GTG".to_vec()]
        );
    }
}

#[test]
fn test_normalize() -> Result<(), LiftOverError> {
    let mut fasta = IndexedReader::from_file(&"testfiles/lift-variant-test/sequence/seq-a.fa")?;

    assert_eq!(
        Variant {
            chromosome: "seq-a".to_string(),
            position: 195,
            reference: b"G".to_vec(),
            alternative: vec![b"T".to_vec()],
        }
        .normalize(&mut fasta)
        .unwrap(),
        Variant {
            chromosome: "seq-a".to_string(),
            position: 195,
            reference: b"G".to_vec(),
            alternative: vec![b"T".to_vec()],
        }
    );

    assert_eq!(
        Variant {
            chromosome: "seq-a".to_string(),
            position: 194,
            reference: b"CGC".to_vec(),
            alternative: vec![b"C".to_vec()],
        }
        .normalize(&mut fasta)
        .unwrap(),
        Variant {
            chromosome: "seq-a".to_string(),
            position: 191,
            reference: b"CCG".to_vec(),
            alternative: vec![b"C".to_vec()],
        }
    );

    assert_eq!(
        Variant {
            chromosome: "seq-a".to_string(),
            position: 11500,
            reference: b"G".to_vec(),
            alternative: vec![b"GTG".to_vec()],
        }
        .normalize(&mut fasta)
        .unwrap(),
        Variant {
            chromosome: "seq-a".to_string(),
            position: 11492,
            reference: b"C".to_vec(),
            alternative: vec![b"CTG".to_vec()],
        }
    );

    assert_eq!(
        Variant {
            chromosome: "seq-a".to_string(),
            position: 11547,
            reference: b"A".to_vec(),
            alternative: vec![b"AAA".to_vec()],
        }
        .normalize(&mut fasta)
        .unwrap(),
        Variant {
            chromosome: "seq-a".to_string(),
            position: 11543,
            reference: b"T".to_vec(),
            alternative: vec![b"TAA".to_vec()],
        }
    );

    assert_eq!(
        Variant {
            chromosome: "seq-a".to_string(),
            position: 11547,
            reference: b"A".to_vec(),
            alternative: vec![b"AAA".to_vec(), b"*".to_vec()],
        }
        .normalize(&mut fasta)
        .unwrap(),
        Variant {
            chromosome: "seq-a".to_string(),
            position: 11543,
            reference: b"T".to_vec(),
            alternative: vec![b"TAA".to_vec(), b"*".to_vec()],
        }
    );

    assert_eq!(
        Variant {
            chromosome: "seq-a".to_string(),
            position: 11547,
            reference: b"AAA".to_vec(),
            alternative: vec![b"A".to_vec()],
        }
        .normalize(&mut fasta)
        .unwrap(),
        Variant {
            chromosome: "seq-a".to_string(),
            position: 11543,
            reference: b"TAA".to_vec(),
            alternative: vec![b"T".to_vec()],
        }
    );

    assert_eq!(
        Variant {
            chromosome: "seq-a".to_string(),
            position: 18040,
            reference: b"GTT".to_vec(),
            alternative: vec![b"GTG".to_vec()],
        }
        .normalize(&mut fasta)
        .unwrap(),
        Variant {
            chromosome: "seq-a".to_string(),
            position: 18042,
            reference: b"T".to_vec(),
            alternative: vec![b"G".to_vec()],
        }
    );

    assert_eq!(
        Variant {
            chromosome: "seq-a".to_string(),
            position: 6655,
            reference: b"ATACCCTAC".to_vec(),
            alternative: vec![b"CTACCCTAA".to_vec()],
        }
        .normalize(&mut fasta)
        .unwrap(),
        Variant {
            chromosome: "seq-a".to_string(),
            position: 6655,
            reference: b"ATACCCTAC".to_vec(),
            alternative: vec![b"CTACCCTAA".to_vec()],
        }
    );

    Ok(())
}

#[test]
fn test_truncate_left_most_nucleotide_if_allele_starts_with_same() {
    let original = Variant {
        chromosome: "hoge".to_string(),
        position: 100,
        reference: b"A".to_vec(),
        alternative: vec![b"C".to_vec()],
    };
    assert_eq!(
        original.truncate_left_most_nucleotide_if_allele_starts_with_same(),
        original
    );

    let original = Variant {
        chromosome: "hoge".to_string(),
        position: 100,
        reference: b"C".to_vec(),
        alternative: vec![b"C".to_vec()],
    };
    assert_eq!(
        original.truncate_left_most_nucleotide_if_allele_starts_with_same(),
        Variant {
            chromosome: "hoge".to_string(),
            position: 101,
            reference: b"".to_vec(),
            alternative: vec![b"".to_vec()],
        }
    );

    let original = Variant {
        chromosome: "hoge".to_string(),
        position: 200,
        reference: b"ATCG".to_vec(),
        alternative: vec![b"ATCA".to_vec()],
    };
    assert_eq!(
        original.truncate_left_most_nucleotide_if_allele_starts_with_same(),
        Variant {
            chromosome: "hoge".to_string(),
            position: 203,
            reference: b"G".to_vec(),
            alternative: vec![b"A".to_vec()],
        }
    );

    let original = Variant {
        chromosome: "hoge".to_string(),
        position: 200,
        reference: b"AG".to_vec(),
        alternative: vec![b"A".to_vec()],
    };
    assert_eq!(
        original.truncate_left_most_nucleotide_if_allele_starts_with_same(),
        Variant {
            chromosome: "hoge".to_string(),
            position: 201,
            reference: b"G".to_vec(),
            alternative: vec![b"".to_vec()],
        }
    );

    let original = Variant {
        chromosome: "hoge".to_string(),
        position: 200,
        reference: b"AG".to_vec(),
        alternative: vec![b"A".to_vec(), b"*".to_vec()],
    };
    assert_eq!(
        original.truncate_left_most_nucleotide_if_allele_starts_with_same(),
        Variant {
            chromosome: "hoge".to_string(),
            position: 201,
            reference: b"G".to_vec(),
            alternative: vec![b"".to_vec(), b"*".to_vec()],
        }
    );
}

#[test]
fn test_truncate_left_most_nucleotide_if_allele_starts_with_same2() {
    assert_eq!(
        Variant {
            chromosome: "hoge".to_string(),
            position: 100,
            reference: b"ATCG".to_vec(),
            alternative: vec![b"CTCA".to_vec()]
        }
        .truncate_left_most_nucleotide_if_allele_starts_with_same(),
        Variant {
            chromosome: "hoge".to_string(),
            position: 100,
            reference: b"ATCG".to_vec(),
            alternative: vec![b"CTCA".to_vec()]
        }
    );

    assert_eq!(
        Variant {
            chromosome: "hoge".to_string(),
            position: 100,
            reference: b"".to_vec(),
            alternative: vec![b"".to_vec()]
        }
        .truncate_left_most_nucleotide_if_allele_starts_with_same(),
        Variant {
            chromosome: "hoge".to_string(),
            position: 100,
            reference: b"".to_vec(),
            alternative: vec![b"".to_vec()]
        }
    );

    assert_eq!(
        Variant {
            chromosome: "hoge".to_string(),
            position: 100,
            reference: b"ATCG".to_vec(),
            alternative: vec![b"ATCA".to_vec()]
        }
        .truncate_left_most_nucleotide_if_allele_starts_with_same(),
        Variant {
            chromosome: "hoge".to_string(),
            position: 103,
            reference: b"G".to_vec(),
            alternative: vec![b"A".to_vec()]
        }
    );

    assert_eq!(
        Variant {
            chromosome: "hoge".to_string(),
            position: 100,
            reference: b"ATCG".to_vec(),
            alternative: vec![b"ATCA".to_vec(), b"*".to_vec()]
        }
        .truncate_left_most_nucleotide_if_allele_starts_with_same(),
        Variant {
            chromosome: "hoge".to_string(),
            position: 103,
            reference: b"G".to_vec(),
            alternative: vec![b"A".to_vec(), b"*".to_vec()]
        }
    );
}
