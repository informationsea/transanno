mod error;

use crate::chain::Strand;
use crate::geneparse::{Feature, FeatureType, Gene, GeneStrand, Transcript};
use crate::poslift::{PositionLiftOver, RegionChangeOp};
pub use error::FeatureLiftError;
use log::debug;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fmt::Write;

impl From<GeneStrand> for Strand {
    fn from(strand: GeneStrand) -> Strand {
        match strand {
            GeneStrand::Forward => Strand::Forward,
            GeneStrand::Reverse => Strand::Reverse,
            GeneStrand::Unknown => unreachable!(),
        }
    }
}

impl From<Strand> for GeneStrand {
    fn from(strand: Strand) -> GeneStrand {
        match strand {
            Strand::Forward => GeneStrand::Forward,
            Strand::Reverse => GeneStrand::Reverse,
        }
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct GeneLiftError<'a, G: Feature, T: Feature, F: Feature> {
    pub error: FeatureLiftError,
    pub gene: &'a Gene<G, T, F>,
    pub failed_transcripts: Vec<TranscriptLiftError<'a, T, F>>,
}

impl<'a, G: Feature, T: Feature, F: Feature> GeneLiftError<'a, G, T, F> {
    pub fn gene_with_failed_reason(&self) -> Gene<G, T, F> {
        let mut new_gene = self.gene.clone();
        new_gene
            .original_record
            .set_attribute("FAIL_REASON", "no_transcripts_lifted");
        new_gene.transcripts = self
            .failed_transcripts
            .iter()
            .map(|x| x.transcript_with_fail_reason())
            .collect();
        new_gene
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct LiftedGene<'a, G: Feature, T: Feature, F: Feature> {
    pub gene: &'a Gene<G, T, F>,
    pub seq_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: GeneStrand,
    pub transcripts: Vec<LiftedTranscript<'a, T, F>>,
    pub failed_transcripts: Vec<TranscriptLiftError<'a, T, F>>,
}

impl<'a, G: Feature, T: Feature, F: Feature> LiftedGene<'a, G, T, F> {
    pub fn apply(&self) -> Gene<G, T, F> {
        let mut new_gene =
            apply_feature_helper(self.gene, &self.seq_id, self.start, self.end, self.strand);
        new_gene.transcripts = self.transcripts.iter().map(|x| x.apply()).collect();
        new_gene
    }

    pub fn gene_with_failed_reason(&self) -> Option<Gene<G, T, F>> {
        if self.failed_transcripts.is_empty() {
            return None;
        }
        let mut new_gene = self.gene.clone();
        new_gene
            .original_record
            .set_attribute("FAIL_REASON", "no_transcripts_lifted");
        new_gene.transcripts = self
            .failed_transcripts
            .iter()
            .map(|x| x.transcript_with_fail_reason())
            .collect();
        Some(new_gene)
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct TranscriptLiftError<'a, T: Feature, F: Feature> {
    pub error: FeatureLiftError,
    pub transcript: &'a Transcript<T, F>,
    pub features: Vec<Vec<LiftedFeature<'a, F>>>,
    pub failed_features: Vec<(&'a F, FeatureLiftError)>,
}

impl<'a, T: Feature, F: Feature> TranscriptLiftError<'a, T, F> {
    pub fn transcript_with_fail_reason(&self) -> Transcript<T, F> {
        let mut new_transcript = self.transcript.clone();
        new_transcript
            .original_record
            .set_attribute("FAIL_REASON", self.error.error_message());
        for one in new_transcript.children.iter_mut() {
            for failed_one in self.failed_features.iter() {
                if one == failed_one.0 {
                    one.set_attribute("FAIL_REASON", failed_one.1.error_message());
                }
            }
            for one_lifted_future in self.features.iter() {
                if one_lifted_future.iter().any(|x| x.feature == one) {
                    one.set_attribute("MAPPED_COUNT", &format!("{}", one_lifted_future.len()));
                }
            }
            if one.attribute("MAPPED_COUNT").is_none() {
                one.set_attribute("MAPPED_COUNT", "0");
            }
        }

        new_transcript
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct LiftedTranscript<'a, T: Feature, F: Feature> {
    pub transcript: &'a Transcript<T, F>,
    pub seq_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: GeneStrand,
    pub features: Vec<LiftedFeature<'a, F>>,
}

impl<'a, T: Feature, F: Feature> LiftedTranscript<'a, T, F> {
    pub fn apply(&self) -> Transcript<T, F> {
        let mut new_transcript = apply_feature_helper(
            self.transcript,
            &self.seq_id,
            self.start,
            self.end,
            self.strand,
        );
        new_transcript.children = self.features.iter().map(|x| x.apply()).collect();
        new_transcript
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct LiftedFeature<'a, F: Feature> {
    pub feature: &'a F,
    pub seq_id: String,
    pub chain_index: usize,
    pub start: u64,
    pub end: u64,
    pub strand: GeneStrand,
    pub changes: Vec<RegionChangeOp>,
}

impl<'a, F: Feature> LiftedFeature<'a, F> {
    pub fn apply(&self) -> F {
        let mut new_feature = apply_feature_helper(
            self.feature,
            &self.seq_id,
            self.start,
            self.end,
            self.strand,
        );

        let ciger = self.changes.iter().fold(String::new(), |mut acc, x| {
            match x {
                RegionChangeOp::Aligned(c) => write!(acc, "M{}", c),
                RegionChangeOp::Deletion(c) => write!(acc, "D{}", c),
                RegionChangeOp::Insertion(c) => write!(acc, "C{}", c),
            }
            .unwrap();
            acc
        });
        new_feature.set_attribute("CIGER", &ciger);

        new_feature
    }
}

fn apply_feature_helper<F: Feature>(
    original_feature: &F,
    new_seq_id: &str,
    new_start: u64,
    new_end: u64,
    new_strand: GeneStrand,
) -> F {
    let mut new_feature = original_feature.clone();
    new_feature.set_attribute("ORIGINAL_CHROM", original_feature.seq_id());
    new_feature.set_attribute("ORIGINAL_START", &format!("{}", original_feature.start()));
    new_feature.set_attribute("ORIGINAL_END", &format!("{}", original_feature.end()));
    new_feature.set_attribute(
        "ORIGINAL_STRAND",
        match original_feature.strand() {
            GeneStrand::Forward => "+",
            GeneStrand::Reverse => "-",
            GeneStrand::Unknown => ".",
        },
    );
    *new_feature.seq_id_mut() = new_seq_id.to_string();
    *new_feature.start_mut() = new_start;
    *new_feature.end_mut() = new_end;
    *new_feature.strand_mut() = new_strand;
    new_feature
}

#[derive(Debug)]
pub struct GeneLiftOver {
    region_lift: PositionLiftOver,
}

impl GeneLiftOver {
    pub fn new(region_lift: PositionLiftOver) -> GeneLiftOver {
        GeneLiftOver { region_lift }
    }

    pub fn lift_single_feature<'a, F: Feature>(
        &self,
        feature: &'a F,
    ) -> Result<Vec<LiftedFeature<'a, F>>, FeatureLiftError> {
        let original_length = feature.len();
        let lower_limit = original_length / 2;
        let upper_limit = original_length * 3 / 2;

        // Lift a region with chain file
        let lifted_region_no_filter = self
            .region_lift
            .lift_region(&feature.seq_id(), feature.range());

        // If length of lifted feature was 50% smaller than original length, mark as failed.
        if lifted_region_no_filter
            .iter()
            .all(|x| lower_limit > x.len())
        {
            if let Some(region) = lifted_region_no_filter
                .iter()
                .find(|x| lower_limit > x.len())
            {
                let e = FeatureLiftError::TooSmallLiftedRegion(region.end - region.start);
                return Err(e);
            }
        }

        // If length of lifted feature was 150% larger than original length, mark as failed.
        if lifted_region_no_filter
            .iter()
            .all(|x| upper_limit < x.len())
        {
            if let Some(region) = lifted_region_no_filter
                .iter()
                .find(|x| upper_limit < x.len())
            {
                let e = FeatureLiftError::TooLargeLiftedRegion(region.len());
                return Err(e);
            }
        }

        let lifted_region = lifted_region_no_filter
            .into_iter()
            .filter(|x| lower_limit <= x.len() && x.len() <= upper_limit)
            .collect::<Vec<_>>();

        // If lift over was failed, mark as failed.
        if lifted_region.is_empty() {
            Err(FeatureLiftError::NoCandidates)
        } else {
            // Multi-map is allowed in this step
            let lifted_features: Vec<_> = lifted_region
                .into_iter()
                .map(|x| {
                    let strand: GeneStrand = x.strand.into();
                    LiftedFeature {
                        feature,
                        seq_id: x.chromosome.name.clone(),
                        chain_index: x.chain_index,
                        start: x.start + 1, // convert to 1-based
                        end: x.end,
                        strand: strand * feature.strand(),
                        changes: x.changes,
                    }
                })
                .collect();
            Ok(lifted_features)
        }
    }

    pub fn lift_gene_feature<'a, G: Feature, T: Feature, F: Feature>(
        &self,
        feature: &'a Gene<G, T, F>,
    ) -> Result<LiftedGene<'a, G, T, F>, GeneLiftError<'a, G, T, F>> {
        // Collect lifted transcripts
        let (success_features, failed_transcripts) = feature.transcripts.iter().fold(
            (Vec::new(), Vec::new()),
            |(mut success, mut failed), x| {
                match self.lift_transcript_feature::<T, F>(x) {
                    Ok(v) => success.push(v),
                    Err(e) => failed.push(e),
                };
                (success, failed)
            },
        );

        if success_features.is_empty() {
            return Err(GeneLiftError {
                error: FeatureLiftError::NoCandidates,
                gene: feature,
                failed_transcripts,
            });
        }

        // Collect all chromosome name from lifted transcripts
        let mut available_chromosomes: Vec<_> = success_features
            .iter()
            .fold(HashMap::new(), |mut acc, x| {
                *acc.entry(x.seq_id.to_string()).or_insert(0) += 1;
                acc
            })
            .into_iter()
            .collect();

        //If two or more chromosomes were found in lifted transcript, pick up most common chromosome.
        available_chromosomes.sort_by(|x, y| match x.1.cmp(&y.1) {
            Ordering::Equal => x.0.cmp(&y.0),
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
        });

        let most_common_chromosome = available_chromosomes.pop().unwrap().0;
        let success_features: Vec<_> = success_features
            .into_iter()
            .filter(|x| x.seq_id == most_common_chromosome)
            .collect();

        // Check strand of all transcripts.
        let mut strand_count: Vec<_> = success_features
            .iter()
            .fold(HashMap::new(), |mut acc, x| {
                *acc.entry(x.strand).or_insert(0) += 1;
                acc
            })
            .into_iter()
            .collect();
        strand_count.sort_by_key(|x| (x.1, x.0));
        let most_common_strand = strand_count.pop().unwrap().0;
        let success_features: Vec<_> = success_features
            .into_iter()
            .filter(|x| x.strand == most_common_strand)
            .collect();

        let new_start = success_features.iter().map(|x| x.start).min().unwrap();
        let new_end = success_features.iter().map(|x| x.end).max().unwrap();

        Ok(LiftedGene {
            gene: feature,
            seq_id: most_common_chromosome,
            start: new_start,
            end: new_end,
            strand: most_common_strand,
            transcripts: success_features,
            failed_transcripts,
        })
    }

    pub fn lift_transcript_feature<'a, T: Feature, F: Feature>(
        &self,
        feature: &'a Transcript<T, F>,
    ) -> Result<LiftedTranscript<'a, T, F>, TranscriptLiftError<'a, T, F>> {
        // Collect lifted exons, CDS, UTR and other features for each transcript.
        let (success_features, failed_features) = feature.children.iter().fold(
            (Vec::new(), Vec::new()),
            |(mut success, mut failed), x| {
                match self.lift_single_feature(x) {
                    Ok(v) => success.push(v),
                    Err(e) => failed.push((x, e)),
                };
                (success, failed)
            },
        );
        if !failed_features.is_empty() || success_features.is_empty() {
            return Err(TranscriptLiftError {
                error: FeatureLiftError::SubFeatureError,
                transcript: feature,
                features: success_features,
                failed_features,
            });
        }

        // Find a chromosome which is included in all lifted features.
        let first_entry_chromosome: HashSet<_> = success_features[0]
            .iter()
            .map(|x| x.seq_id.clone())
            .collect();
        let common_chromosome =
            success_features
                .iter()
                .fold(first_entry_chromosome, |mut acc, x| {
                    let this_entry_chromosome: HashSet<_> = x.iter().map(|x| &x.seq_id).collect();
                    acc.retain(|x| this_entry_chromosome.contains(x));
                    acc
                });
        // If two or more chromosomes were found, mark this transcript as failed.
        if common_chromosome.len() >= 2 {
            return Err(TranscriptLiftError {
                error: FeatureLiftError::MultiMap,
                transcript: feature,
                features: success_features,
                failed_features,
            });
        }
        // If no chromosome was found, mark this transcript as failed.
        if common_chromosome.is_empty() {
            return Err(TranscriptLiftError {
                error: FeatureLiftError::NoCommonChromosome,
                transcript: feature,
                features: success_features,
                failed_features,
            });
        }
        // If a part of exons were lifted to multiple chromosome, remove lifted region which is not lifted to selected chromosome from candidates.
        let mut success_features: Vec<_> = success_features
            .into_iter()
            .map(|x| {
                x.into_iter()
                    .filter(|y| common_chromosome.contains(&y.seq_id))
                    .collect::<Vec<_>>()
            })
            .collect();
        if let Some(_zero) = success_features.iter().find(|x| x.is_empty()) {
            unreachable!();
        }
        // Find common chain index
        let mut common_chain_index: HashSet<_> =
            success_features[0].iter().map(|x| x.chain_index).collect();
        for one_success_futures in &success_features {
            common_chain_index.retain(|x| one_success_futures.iter().any(|y| y.chain_index == *x));
        }

        // Remove all features which are not included in common chain index, if common chain index is found.
        if common_chain_index.len() == 1 && success_features.iter().any(|x| x.len() > 1) {
            success_features = success_features
                .into_iter()
                .map(|x| {
                    x.into_iter()
                        .filter(|y| common_chain_index.contains(&y.chain_index))
                        .collect::<Vec<_>>()
                })
                .collect();
        }

        // Check multi-mapped features
        if let Some(_multi) = success_features.iter().find(|x| x.len() >= 2) {
            return Err(TranscriptLiftError {
                error: FeatureLiftError::MultiMapSubFeature,
                transcript: feature,
                features: success_features,
                failed_features,
            });
        }
        let success_features: Vec<_> = success_features
            .into_iter()
            .map(|mut x| x.pop().unwrap())
            .collect();

        // Check all features have same strand
        if success_features
            .iter()
            .any(|x| success_features[0].strand != x.strand)
        {
            debug!("wrong strand: {:?}", feature);
            return Err(TranscriptLiftError {
                error: FeatureLiftError::WrongStrand,
                transcript: feature,
                features: success_features.into_iter().map(|x| vec![x]).collect(),
                failed_features,
            });
        }

        // Check the order of exons
        let mut original_exons: Vec<_> = feature
            .children
            .iter()
            .filter(|x| x.feature_type() == FeatureType::Exon)
            .collect();
        original_exons.sort_by_key(|x| x.start());

        let mut lifted_exons: Vec<_> = success_features
            .iter()
            .filter(|x| x.feature.feature_type() == FeatureType::Exon)
            .collect();
        lifted_exons.sort_by_key(|x| x.start);

        if original_exons[0].strand() != lifted_exons[0].strand {
            lifted_exons.reverse();
        }

        let lifted_exon_id_list: Vec<_> = lifted_exons.iter().map(|x| x.feature).collect();
        if original_exons != lifted_exon_id_list {
            debug!("wrong exon order: {:?}", feature);
            return Err(TranscriptLiftError {
                error: FeatureLiftError::WrongExonOrder,
                transcript: feature,
                features: success_features.into_iter().map(|x| vec![x]).collect(),
                failed_features,
            });
        }

        // Define new region from features.
        let new_start = success_features.iter().map(|x| x.start).min().unwrap();
        let new_end = success_features.iter().map(|x| x.end).max().unwrap();

        let upper_limit = (feature.end() - feature.start()) * 3 / 2;
        let lower_limit = (feature.end() - feature.start()) / 2;
        // If length of transcript was 150% larger than original length, mark as failed.
        if upper_limit < new_end - new_start {
            return Err(TranscriptLiftError {
                error: FeatureLiftError::TooLargeLiftedRegion(new_end - new_start),
                transcript: feature,
                features: success_features.into_iter().map(|x| vec![x]).collect(),
                failed_features,
            });
        }
        // If length of transcript was 50% smaller than original length, mark as failed.
        if lower_limit > new_end - new_start {
            return Err(TranscriptLiftError {
                error: FeatureLiftError::TooSmallLiftedRegion(new_end - new_start),
                transcript: feature,
                features: success_features.into_iter().map(|x| vec![x]).collect(),
                failed_features,
            });
        }

        Ok(LiftedTranscript {
            seq_id: success_features[0].seq_id.clone(),
            strand: success_features[0].strand,
            features: success_features,
            start: new_start,
            end: new_end,
            transcript: feature,
        })
    }
}

#[cfg(test)]
mod test;
