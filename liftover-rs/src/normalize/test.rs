use super::*;
use bio::io::fasta::IndexedReader;
use serde::Deserialize;

#[allow(clippy::unreadable_literal)]
#[test]
fn test_extend_to_left_if_empty_allele_exists() -> anyhow::Result<()> {
    let mut fasta = IndexedReader::from_file(&"testfiles/genomes/GRCh38/GRCh38.chr22.genome.fa")?;
    {
        let chromosome = "chr22";
        let mut position = 17513945;
        let mut alleles = vec![b"".to_vec(), b"G".to_vec(), b"GC".to_vec()];
        assert!(extend_to_left_if_empty_allele_exists(
            chromosome,
            &mut fasta,
            &mut position,
            &mut alleles
        )?);
        assert_eq!(position, 17513944);
        assert_eq!(
            alleles,
            vec![b"A".to_vec(), b"AG".to_vec(), b"AGC".to_vec()]
        );
    }

    {
        let chromosome = "chr22";
        let mut position = 17513945;
        let mut alleles = vec![b"".to_vec(), b"G".to_vec(), b"GC".to_vec(), b"*".to_vec()];
        assert!(extend_to_left_if_empty_allele_exists(
            chromosome,
            &mut fasta,
            &mut position,
            &mut alleles
        )?);
        assert_eq!(position, 17513944);
        assert_eq!(
            alleles,
            vec![
                b"A".to_vec(),
                b"AG".to_vec(),
                b"AGC".to_vec(),
                b"*".to_vec()
            ]
        );
    }

    {
        let chromosome = "chr22";
        let mut position = 17513944;
        let mut alleles = vec![b"A".to_vec(), b"AG".to_vec(), b"AGC".to_vec()];
        assert!(!extend_to_left_if_empty_allele_exists(
            chromosome,
            &mut fasta,
            &mut position,
            &mut alleles
        )?);
        assert_eq!(position, 17513944);
        assert_eq!(
            alleles,
            vec![b"A".to_vec(), b"AG".to_vec(), b"AGC".to_vec()],
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

#[derive(Deserialize, Debug, PartialEq, Eq)]
struct NormalizeTestSet {
    chromosome: String,
    original_pos: u64,
    original_reference: String,
    original_alternative: String,
    expected_pos: u64,
    expected_reference: String,
    expected_alternative: String,
}

#[test]
fn test_normalize() -> anyhow::Result<()> {
    let mut fasta = IndexedReader::from_file(&"testfiles/genomes/GRCh38/GRCh38.chr22.genome.fa")?;
    let testset = include_bytes!("normalize-testset.csv"); // 1-based
    let mut rdr = csv::Reader::from_reader(&testset[..]);
    for result in rdr.deserialize() {
        let record: NormalizeTestSet = result?;
        //println!("record {:?}", record);
        let normalized = Variant::new(
            &record.chromosome,
            record.original_pos - 1, // convert to 0-based
            &record.original_reference.as_bytes(),
            &record
                .original_alternative
                .split(',')
                .map(|x| x.as_bytes())
                .collect::<Vec<_>>(),
        )
        .normalize(&mut fasta)?;
        //println!("normalized: {:?}", normalized);
        assert_eq!(normalized.position, record.expected_pos - 1); // convert to 0-based
    }

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
