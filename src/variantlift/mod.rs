use crate::chain::*;
use crate::poslift::*;
use crate::{reverse_complement, GenomeSequence, LiftOverError, Variant};
use std::io::prelude::*;
use std::str;

pub mod error;

pub type VariantLiftOverResult = Result<LiftedVariant, error::VariantLiftOverError>;

#[derive(Debug, Hash, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct LiftedVariant {
    pub chromosome: String,
    pub position: u64,
    pub strand: Strand,
    pub original_reference: Vec<u8>,
    pub reference: Vec<u8>,
    pub alternative: Vec<Vec<u8>>,
    pub reference_changed: bool,
}

impl<'a> From<LiftedVariant> for Variant {
    fn from(v: LiftedVariant) -> Variant {
        Variant {
            chromosome: v.chromosome,
            position: v.position,
            reference: v.reference,
            alternative: v.alternative,
        }
    }
}

pub trait LiftVariant {
    fn lift_variant(
        &mut self,
        variant: &Variant,
        acceptable_deletion: u64,
        acceptable_insertion: u64,
    ) -> Result<Vec<VariantLiftOverResult>, LiftOverError>;
}

pub struct VariantLiftOver<G: GenomeSequence, L: LiftPosition> {
    reference_sequence: G,
    query_sequence: G,
    lift_position: L,
}

impl<G: GenomeSequence> VariantLiftOver<G, PositionLiftOver> {
    pub fn load<R: Read>(
        chain_file: R,
        reference_sequence: G,
        query_sequence: G,
    ) -> Result<Self, LiftOverError> {
        let lift_position = PositionLiftOver::load(chain_file)?;
        Ok(VariantLiftOver {
            reference_sequence,
            query_sequence,
            lift_position,
        })
    }

    pub fn new(chain_file: ChainFile, reference_sequence: G, query_sequence: G) -> Self {
        let lift_position = PositionLiftOver::new(chain_file);
        VariantLiftOver {
            reference_sequence,
            query_sequence,
            lift_position,
        }
    }
}

impl<G: GenomeSequence, L: LiftPosition> LiftPosition for VariantLiftOver<G, L> {
    fn chain_list(&self) -> &[Chain] {
        self.lift_position.chain_list()
    }
    fn reference_chromosomes(&self) -> &[Chromosome] {
        self.lift_position.reference_chromosomes()
    }

    fn query_chromosomes(&self) -> &[Chromosome] {
        self.lift_position.query_chromosomes()
    }

    fn reference_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.lift_position.reference_chromosome_by_name(name)
    }

    fn query_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.lift_position.query_chromosome_by_name(name)
    }

    fn lift_position<'a>(&'a self, chromosome: &str, position: u64) -> Vec<LiftOverResult> {
        self.lift_position.lift_position(chromosome, position)
    }
    fn search_target<'a>(
        &'a self,
        chromosome: &str,
        start: u64,
        end: u64,
    ) -> Vec<Vec<&'a TargetRegion>> {
        self.lift_position.search_target(chromosome, start, end)
    }

    fn lift_region(&self, chromosome: &str, start: u64, end: u64) -> Vec<LiftRegionResult> {
        self.lift_position.lift_region(chromosome, start, end)
    }
}

impl<G: GenomeSequence, L: LiftPosition> LiftVariant for VariantLiftOver<G, L> {
    // Input variant should be left aligned
    fn lift_variant(
        &mut self,
        variant: &Variant,
        acceptable_deletion: u64,
        acceptable_insertion: u64,
    ) -> Result<Vec<VariantLiftOverResult>, LiftOverError> {
        // check reference sequence name
        if self
            .reference_chromosome_by_name(&variant.chromosome)
            .is_none()
        {
            return Ok(vec![Err(error::VariantLiftOverError::UnknownSequenceName(
                variant.chromosome.to_string(),
            ))]);
        }

        let start = variant.position;
        let end = variant.position + variant.reference.len() as u64;
        let largest_alternative = variant
            .alternative
            .iter()
            .map(|x| x.len())
            .max()
            .unwrap_or(0);

        let expected_ref = self
            .reference_sequence
            .get_sequence(&variant.chromosome, start, end)?;
        if expected_ref != variant.reference {
            return Ok(vec![Err(
                error::VariantLiftOverError::ReferenceSequenceIsNotMatch,
            )]);
        }

        let target_region_lists =
            if variant.reference.len() == 1 && variant.alternative.iter().all(|x| x.len() == 1) {
                // do not search neighborhood variants for SNPs
                self.lift_position
                    .search_target(&variant.chromosome, start, end)
            } else {
                self.lift_position
                    .search_target(&variant.chromosome, start - 1, end + 1)
            };

        let mut result = Vec::new();

        for one_target_list in target_region_lists {
            let consider_start = one_target_list
                .iter()
                .filter(|x| x.is_in_gap)
                .map(|x| x.reference_start)
                .chain(vec![start].into_iter())
                .min()
                .unwrap();
            let consider_end = one_target_list
                .iter()
                .filter(|x| x.is_in_gap)
                .map(|x| x.reference_end + 1)
                .chain(vec![end].into_iter())
                .max()
                .unwrap();
            // println!(
            //     "consider: {}:{}-{}",
            //     variant.chromosome, consider_start, consider_end
            // );
            let one_region = if let Some(v) = self
                .lift_position
                .lift_region(&variant.chromosome, consider_start, consider_end)
                .into_iter()
                .find(|x| x.chain_index == one_target_list[0].chain_index)
            {
                v
            } else {
                // corresponding region was not found
                continue;
            };

            let is_acceptable_deletion = |x: &TargetRegion| {
                if x.is_in_gap {
                    x.reference_len() <= acceptable_deletion + variant.reference.len() as u64
                } else {
                    true
                }
            };
            let is_acceptable_insertion = |x: &TargetRegion| {
                if x.is_in_gap {
                    x.query_len() <= acceptable_insertion + largest_alternative as u64
                } else {
                    true
                }
            };

            if one_target_list.iter().any(|x| !is_acceptable_deletion(x)) {
                result.push(Err(
                    error::VariantLiftOverError::UnacceptableLargeDeletion {
                        chromosome: one_region.chromosome.name.to_string(),
                        start: one_region.start,
                        end: one_region.end,
                    },
                ));
                continue;
            }

            if one_target_list.iter().any(|x| !is_acceptable_insertion(x)) {
                result.push(Err(
                    error::VariantLiftOverError::UnacceptableLargeInsertion {
                        chromosome: one_region.chromosome.name.to_string(),
                        start: one_region.start,
                        end: one_region.end,
                    },
                ));
                continue;
            }

            let extended_start_seq =
                self.reference_sequence
                    .get_sequence(&variant.chromosome, consider_start, start)?;
            let extended_end_seq =
                self.reference_sequence
                    .get_sequence(&variant.chromosome, end, consider_end)?;

            let query_chromosome = one_region.chromosome;
            let query_sequence = self.query_sequence.get_sequence(
                &one_region.chromosome.name,
                one_region.start,
                one_region.end,
            )?;

            let reference_seq_original: Vec<u8> = extended_start_seq
                .iter()
                .chain(variant.reference.iter())
                .chain(extended_end_seq.iter())
                .cloned()
                .collect();
            let alternate_seq_original: Vec<Vec<u8>> = variant
                .alternative
                .iter()
                .map(|x| {
                    if x == b"*" {
                        b"*".to_vec()
                    } else {
                        extended_start_seq
                            .iter()
                            .chain(x.iter())
                            .chain(extended_end_seq.iter())
                            .cloned()
                            .collect()
                    }
                })
                .collect();

            let reference_seq = match one_region.strand {
                Strand::Forward => reference_seq_original,
                Strand::Reverse => reverse_complement(&reference_seq_original),
            };
            let alternate_seq = match one_region.strand {
                Strand::Forward => alternate_seq_original,
                Strand::Reverse => alternate_seq_original
                    .iter()
                    .map(|x| reverse_complement(&x))
                    .collect(),
            };

            let mut variant = Variant {
                chromosome: query_chromosome.name.to_string(),
                position: one_region.start,
                reference: query_sequence,
                alternative: vec![reference_seq]
                    .into_iter()
                    .chain(alternate_seq.into_iter())
                    .collect(),
            }
            .normalize(&mut self.query_sequence)?;

            let original_reference = variant.alternative.remove(0);

            result.push(Ok(LiftedVariant {
                chromosome: variant.chromosome,
                position: variant.position,
                strand: one_region.strand,
                reference_changed: original_reference != variant.reference,
                original_reference,
                reference: variant.reference,
                alternative: variant.alternative,
            }));
        }

        Ok(result)
    }
}

#[cfg(test)]
mod test;
