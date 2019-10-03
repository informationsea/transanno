use crate::chain::*;
use crate::LiftOverError;
use bio::data_structures::interval_tree::IntervalTree;
use log::trace;
use regex::Regex;

use std::collections::HashMap;
use std::io::Read;
use std::u64;

lazy_static! {
    static ref SPACE: Regex = Regex::new("\\s").unwrap();
}

#[derive(Hash, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Copy)]
pub enum RegionChangeOp {
    Aligned(u64),
    Insertion(u64),
    Deletion(u64),
}

#[derive(Hash, Debug, Clone, PartialEq, Eq)]
pub struct LiftRegionResult<'a> {
    pub chromosome: &'a Chromosome,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub changes: Vec<RegionChangeOp>,
    pub chain_index: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetRegion {
    pub reference_chromosome_index: usize,
    pub query_chromosome_index: usize,
    pub chain_index: usize,
    pub reference_start: u64,
    pub reference_end: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: Strand,
    pub is_in_gap: bool,
}

impl TargetRegion {
    pub fn contains_reference_position(&self, position: u64) -> bool {
        self.reference_start <= position && position < self.reference_end
    }

    pub fn contains_query_position(&self, position: u64) -> bool {
        self.query_start <= position && position < self.query_end
    }

    pub fn is_indel(&self) -> bool {
        (self.reference_end - self.reference_start) != (self.query_end - self.query_start)
    }

    pub fn reference_len(&self) -> u64 {
        self.reference_end - self.reference_start
    }

    pub fn query_len(&self) -> u64 {
        self.query_end - self.query_start
    }
}

#[derive(Debug, Hash, Clone, PartialEq)]
pub struct LiftOverResult<'a> {
    pub chromosome: &'a Chromosome,
    pub position: u64,
    pub chain_index: usize,
    pub strand: Strand,
}

pub trait LiftPosition {
    fn chain_list(&self) -> &[Chain];
    fn reference_chromosomes(&self) -> &[Chromosome];
    fn query_chromosomes(&self) -> &[Chromosome];
    fn reference_chromosome_by_name(&self, name: &str) -> Option<&Chromosome>;
    fn query_chromosome_by_name(&self, name: &str) -> Option<&Chromosome>;
    fn lift_position(&self, chromosome: &str, position: u64) -> Vec<LiftOverResult>;
    fn search_target<'a>(
        &'a self,
        chromosome: &str,
        start: u64,
        end: u64,
    ) -> Vec<Vec<&'a TargetRegion>>;
    fn lift_region<'a>(
        &'a self,
        chromosome: &str,
        start: u64,
        end: u64,
    ) -> Vec<LiftRegionResult<'a>>;
}

#[derive(Debug)]
pub struct PositionLiftOver {
    reference_interval: HashMap<String, IntervalTree<u64, TargetRegion>>,
    query_interval: HashMap<String, IntervalTree<u64, TargetRegion>>,
    chain_file: ChainFile,
}

impl LiftPosition for PositionLiftOver {
    fn chain_list(&self) -> &[Chain] {
        &self.chain_file.chain_list
    }

    fn reference_chromosomes(&self) -> &[Chromosome] {
        &self.chain_file.reference_chromosomes
    }

    fn query_chromosomes(&self) -> &[Chromosome] {
        &self.chain_file.query_chromosomes
    }

    fn reference_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.chain_file.reference_chromosome_by_name(name)
    }

    fn query_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.chain_file.query_chromosome_by_name(name)
    }

    // LiftOver a genomic coordinate. A position should be zero-based.
    fn lift_position(&self, chromosome: &str, position: u64) -> Vec<LiftOverResult> {
        let mut results = Vec::new();
        for loop_one in self.search_target(chromosome, position, position + 1) {
            for data in loop_one {
                // skip indel
                if (data.reference_end - data.reference_start)
                    != (data.query_end - data.query_start)
                {
                    continue;
                }

                if data.is_in_gap {
                    continue;
                }

                results.push(LiftOverResult {
                    chromosome: self
                        .query_chromosomes()
                        .get(data.query_chromosome_index)
                        .unwrap(),
                    position: match data.strand {
                        Strand::Forward => data.query_start + (position - data.reference_start),
                        Strand::Reverse => data.query_start + (data.reference_end - position) - 1,
                    },
                    chain_index: data.chain_index,
                    strand: data.strand,
                });
            }
        }
        results
    }

    fn search_target<'a>(
        &'a self,
        chromosome: &str,
        start: u64,
        end: u64,
    ) -> Vec<Vec<&'a TargetRegion>> {
        let mut groups = HashMap::new();
        if let Some(interval_tree) = self.get_reference_interval(chromosome) {
            for one in interval_tree.find(start..end) {
                let data: &TargetRegion = one.data();
                if data.reference_end == data.reference_start && start == data.reference_start {
                    continue; // skip insertion at head of region
                }
                groups
                    .entry(data.chain_index)
                    .or_insert_with(Vec::new)
                    .push(data);
            }
        }

        for one_group in groups.values_mut() {
            one_group.sort_by_key(|x| match x.strand {
                Strand::Forward => (
                    x.reference_start,
                    x.reference_end,
                    x.query_start,
                    x.query_end,
                ),
                Strand::Reverse => (
                    x.reference_start,
                    x.reference_end,
                    u64::MAX - x.query_start,
                    u64::MAX - x.query_end,
                ),
            });
        }

        groups.into_iter().map(|x| x.1).collect()
    }

    fn lift_region(&self, chromosome: &str, start: u64, end: u64) -> Vec<LiftRegionResult> {
        let mut result = Vec::new();
        for one_chain in self.search_target(chromosome, start, end) {
            // println!("lift region: {:?}", one_chain);
            if one_chain.len() == 1 {
                let one_target = one_chain[0];
                // println!(
                //     "check one target {} {} {}",
                //     one_target.contains_reference_position(start),
                //     one_target.contains_reference_position(end),
                //     !one_target.is_indel()
                // );
                if one_target.reference_start <= start
                    && end <= one_target.reference_end
                    && !one_target.is_in_gap
                {
                    let start_offset = start - one_target.reference_start;
                    let end_offset = one_target.reference_end - end;

                    result.push(LiftRegionResult {
                        chromosome: &self.query_chromosomes()[one_chain[0].query_chromosome_index],
                        start: match one_target.strand {
                            Strand::Forward => one_target.query_start + start_offset,
                            Strand::Reverse => one_target.query_start + end_offset,
                        },
                        end: match one_target.strand {
                            Strand::Forward => one_target.query_end - end_offset,
                            Strand::Reverse => one_target.query_end - start_offset,
                        },
                        strand: one_target.strand,
                        changes: vec![RegionChangeOp::Aligned(end - start)],
                        chain_index: one_target.chain_index,
                    });
                } else {
                    //println!("skipped one target");
                    // out of chain
                    continue;
                }
            } else {
                // multiple target region
                let mut changes = Vec::new();

                // first target
                if one_chain[0].reference_start > start {
                    //println!("skipped first target");
                    continue; // out of chain
                }

                let first_offset = start - one_chain[0].reference_start;
                let query_start = if !one_chain[0].is_in_gap {
                    // aligned region
                    changes.push(RegionChangeOp::Aligned(
                        one_chain[0].reference_end - one_chain[0].reference_start - first_offset,
                    ));

                    match one_chain[0].strand {
                        Strand::Forward => one_chain[0].query_start + first_offset,
                        Strand::Reverse => one_chain[0].query_end - first_offset,
                    }
                } else {
                    // INDEL
                    changes.push(RegionChangeOp::Deletion(
                        one_chain[0].reference_end - one_chain[0].reference_start - first_offset,
                    ));

                    if one_chain[1].is_in_gap {
                        // INDEL should not be repeated
                        unreachable!()
                    }

                    match one_chain[0].strand {
                        Strand::Forward => one_chain[0].query_end,
                        Strand::Reverse => one_chain[0].query_start,
                    }
                };

                for one_target in &one_chain[1..(one_chain.len() - 1)] {
                    if !one_target.is_in_gap {
                        // aligned region
                        changes.push(RegionChangeOp::Aligned(
                            one_target.reference_end - one_target.reference_start,
                        ));
                    } else {
                        let insert_length = one_target.query_end - one_target.query_start;
                        let delete_length = one_target.reference_end - one_target.reference_start;
                        if insert_length > 0 {
                            changes.push(RegionChangeOp::Insertion(insert_length));
                        }
                        if delete_length > 0 {
                            changes.push(RegionChangeOp::Deletion(delete_length));
                        }
                    }
                }

                // last target
                let last_target = one_chain.last().unwrap();
                if last_target.reference_end < end {
                    // out of chain
                    //println!("skipped last target");
                    continue;
                }
                let last_offset = last_target.reference_end - end;
                let query_end = if !last_target.is_indel() {
                    // aligned region
                    changes.push(RegionChangeOp::Aligned(
                        last_target.reference_end - last_target.reference_start - last_offset,
                    ));

                    match one_chain[0].strand {
                        Strand::Forward => last_target.query_end - last_offset,
                        Strand::Reverse => last_target.query_start + last_offset,
                    }
                } else {
                    // INDEL
                    changes.push(RegionChangeOp::Deletion(
                        last_target.reference_end - last_target.reference_start - last_offset,
                    ));

                    match one_chain[0].strand {
                        Strand::Forward => last_target.query_start,
                        Strand::Reverse => last_target.query_end,
                    }
                };

                result.push(LiftRegionResult {
                    chromosome: &self.query_chromosomes()[one_chain[0].query_chromosome_index],
                    start: match one_chain[0].strand {
                        Strand::Forward => query_start,
                        Strand::Reverse => query_end,
                    },
                    end: match one_chain[0].strand {
                        Strand::Forward => query_end,
                        Strand::Reverse => query_start,
                    },
                    strand: one_chain[0].strand,
                    changes,
                    chain_index: last_target.chain_index,
                });
            }
        }
        result
    }
}

fn chain_register_interval_tree(
    chain: &Chain,
    reference_interval: &mut IntervalTree<u64, TargetRegion>,
    query_interval: &mut IntervalTree<u64, TargetRegion>,
    reference_chromosome_index: usize,
    query_chromosome_index: usize,
    chain_index: usize,
) {
    assert_eq!(chain.reference_strand, Strand::Forward);

    let mut reference_current = chain.reference_start;
    let mut query_current = chain.query_start;
    for one in chain.chain_interval.iter() {
        // insert un-gapped region
        let reference_next = reference_current + one.size;
        let query_next = query_current + one.size;

        register_one_interval(
            chain,
            reference_interval,
            query_interval,
            reference_chromosome_index,
            query_chromosome_index,
            chain_index,
            reference_current,
            reference_next,
            query_current,
            query_next,
            false,
        );

        reference_current = reference_next;
        query_current = query_next;

        // insert gapped region
        if let Some(diff_ref) = one.difference_reference {
            if let Some(diff_query) = one.difference_query {
                let reference_next = reference_current + diff_ref;
                let query_next = query_current + diff_query;

                register_one_interval(
                    chain,
                    reference_interval,
                    query_interval,
                    reference_chromosome_index,
                    query_chromosome_index,
                    chain_index,
                    reference_current,
                    reference_next,
                    query_current,
                    query_next,
                    true,
                );

                reference_current = reference_next;
                query_current = query_next;
            } else {
                panic!("Invalid format")
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn register_one_interval(
    chain: &Chain,
    reference_interval: &mut IntervalTree<u64, TargetRegion>,
    query_interval: &mut IntervalTree<u64, TargetRegion>,
    reference_chromosome_index: usize,
    query_chromosome_index: usize,
    chain_index: usize,
    reference_current: u64,
    reference_next: u64,
    query_current: u64,
    query_next: u64,
    is_in_gap: bool,
) {
    match chain.query_strand {
        Strand::Forward => {
            let target_region = TargetRegion {
                reference_chromosome_index,
                query_chromosome_index,
                chain_index,
                reference_start: reference_current,
                reference_end: reference_next,
                query_start: query_current,
                query_end: query_next,
                strand: Strand::Forward,
                is_in_gap,
            };

            let register_reference_next = if reference_current == reference_next {
                reference_next + 1
            } else {
                reference_next
            };
            let register_query_next = if query_current == query_next {
                query_next + 1
            } else {
                query_next
            };

            reference_interval.insert(
                reference_current..register_reference_next,
                target_region.clone(),
            );
            query_interval.insert(query_current..register_query_next, target_region);
        }
        Strand::Reverse => {
            let target_region = TargetRegion {
                reference_chromosome_index,
                query_chromosome_index,
                chain_index,
                reference_start: reference_current,
                reference_end: reference_next,
                query_start: chain.query_chromosome.length - query_next,
                query_end: chain.query_chromosome.length - query_current,
                strand: Strand::Reverse,
                is_in_gap,
            };

            let register_reference_next = if reference_current == reference_next {
                reference_next + 1
            } else {
                reference_next
            };
            let register_query_next = if query_current == query_next {
                query_next + 1
            } else {
                query_next
            };

            reference_interval.insert(
                reference_current..register_reference_next,
                target_region.clone(),
            );
            query_interval.insert(
                (chain.query_chromosome.length - register_query_next)
                    ..(chain.query_chromosome.length - query_current),
                target_region,
            );
        }
    }
}

impl PositionLiftOver {
    pub fn new(chain_file: ChainFile) -> PositionLiftOver {
        let mut reference_interval = HashMap::new();
        let mut query_interval = HashMap::new();

        for (i, one_chain) in chain_file.chain_list.iter().enumerate() {
            trace!(
                "Register chain: {}  ({} -> {})",
                one_chain.chain_id,
                one_chain.reference_chromosome.name,
                one_chain.query_chromosome.name
            );

            let one_reference_interval = reference_interval
                .entry(one_chain.reference_chromosome.name.clone())
                .or_insert_with(IntervalTree::new);
            let one_query_interval = query_interval
                .entry(one_chain.query_chromosome.name.clone())
                .or_insert_with(IntervalTree::new);

            chain_register_interval_tree(
                one_chain,
                one_reference_interval,
                one_query_interval,
                chain_file.reference_chromosome_name_to_index[&one_chain.reference_chromosome.name],
                chain_file.query_chromosome_name_to_index[&one_chain.query_chromosome.name],
                i,
            );
        }

        PositionLiftOver {
            reference_interval,
            query_interval,
            chain_file,
        }
    }

    pub fn load<R: Read>(chain_reader: R) -> Result<PositionLiftOver, LiftOverError> {
        let chain_file = ChainFile::load(chain_reader)?;
        Ok(PositionLiftOver::new(chain_file))
    }

    fn get_reference_interval(&self, chromosome: &str) -> Option<&IntervalTree<u64, TargetRegion>> {
        if let Some(val) = self.reference_interval.get(chromosome) {
            Some(val)
        } else if let Some(val) = self.reference_interval.get(&format!("chr{}", chromosome)) {
            Some(val)
        } else if chromosome.len() >= 4 && chromosome.starts_with("chr") {
            if let Some(val) = self.reference_interval.get(&chromosome[3..]) {
                Some(val)
            } else {
                None
            }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test;
