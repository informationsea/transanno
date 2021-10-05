use super::*;
use serde::Deserialize;
use std::collections::HashSet;

#[derive(Clone, Debug, Deserialize, PartialEq, Eq)]
struct LiftVariantExpected {
    original_chrom: String,
    original_pos: u64,
    original_ref: String,
    original_alt: String,
    mapped_chrom: String,
    mapped_pos: String,
    mapped_original_ref: String,
    mapped_ref: String,
    mapped_alt: String,
    mapped_ref_changed: String,
    mapped_strand: String,
    comment: String,
}

impl LiftVariantExpected {
    pub fn original_variant(&self) -> Variant {
        Variant::new(
            &self.original_chrom,
            self.original_pos - 1,
            self.original_ref.as_bytes(),
            &self
                .original_alt
                .split(',')
                .map(|x| x.as_bytes())
                .collect::<Vec<_>>(),
        )
    }

    pub fn mapped_variant(&self) -> HashSet<LiftedVariant> {
        self.mapped_chrom
            .split('|')
            .zip(self.mapped_pos.split('|'))
            .zip(self.mapped_ref.split('|'))
            .zip(self.mapped_alt.split('|'))
            .zip(self.mapped_original_ref.split('|'))
            .zip(self.mapped_ref_changed.split('|'))
            .zip(self.mapped_strand.split('|'))
            .map(
                |(
                    (
                        ((((one_chrom, one_pos), one_ref), one_alt), one_mapped_original_ref),
                        one_mapped_ref_changed,
                    ),
                    one_mapped_strand,
                )| {
                    LiftedVariant {
                        chromosome: one_chrom.to_string(),
                        position: one_pos.parse::<u64>().unwrap() - 1,
                        original_reference: one_mapped_original_ref.as_bytes().to_vec(),
                        reference: one_ref.as_bytes().to_vec(),
                        alternative: one_alt
                            .split(',')
                            .map(|x| x.as_bytes().to_vec())
                            .collect::<Vec<_>>(),
                        reference_changed: one_mapped_ref_changed == "TRUE",
                        strand: if one_mapped_strand == "+" {
                            Strand::Forward
                        } else {
                            Strand::Reverse
                        },
                    }
                },
            )
            .collect()
    }
}

#[test]
fn test_lift_variant() -> anyhow::Result<()> {
    let mut reference_fasta = bio::io::fasta::IndexedReader::from_file(
        &"testfiles/genomes/GRCh38/GRCh38.chr22.genome.fa",
    )?;
    let mut query_fasta = bio::io::fasta::IndexedReader::from_file(
        &"testfiles/genomes/GRCh37/GRCh37.chr22.genome.fa",
    )?;

    let chain_file = ChainFile::load(
        &include_bytes!("../../testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.chain")[..],
    )?
    .left_align(&mut reference_fasta, &mut query_fasta)?;

    let mut variant_liftover = VariantLiftOver::new(chain_file, reference_fasta, query_fasta);

    let mut expected_reader = csv::Reader::from_reader(&include_bytes!("lift-test.csv")[..]);
    for expected in expected_reader.deserialize() {
        let expected: LiftVariantExpected = expected?;
        if expected.original_chrom.starts_with('#') {
            continue;
        }
        let original = expected.original_variant();
        let expected_variants = expected.mapped_variant();
        let result = variant_liftover.lift_variant(&original, 10, 10)?;
        //println!("expected: {:?}   result: {:?}", expected_variants, result);
        assert_eq!(
            expected_variants,
            result
                .iter()
                .map(|x| x.as_ref().unwrap())
                .cloned()
                .collect()
        );
    }

    Ok(())
}
