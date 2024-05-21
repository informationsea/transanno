use crate::chain::*;
use crate::LiftOverError;
use bio::data_structures::interval_tree::IntervalTree;
use log::trace;

use std::collections::HashMap;
use std::io::Read;
use std::ops::Range;
use std::u64;

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

impl<'a> LiftRegionResult<'a> {
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }

    pub fn range(&self) -> std::ops::Range<u64> {
        self.start..self.end
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetRegion {
    pub original_chromosome_index: usize,
    pub new_chromosome_index: usize,
    pub chain_index: usize,
    pub original_start: u64,
    pub original_end: u64,
    pub new_start: u64,
    pub new_end: u64,
    pub strand: Strand,
    pub is_in_gap: bool,
}

impl TargetRegion {
    pub fn contains_original_position(&self, position: u64) -> bool {
        self.original_start <= position && position < self.original_end
    }

    pub fn contains_new_position(&self, position: u64) -> bool {
        self.new_start <= position && position < self.new_end
    }

    pub fn is_indel(&self) -> bool {
        (self.original_end - self.original_start) != (self.new_end - self.new_start)
    }

    pub fn original_len(&self) -> u64 {
        self.original_end - self.original_start
    }

    pub fn new_len(&self) -> u64 {
        self.new_end - self.new_start
    }
}

#[derive(Debug, Hash, Clone, PartialEq)]
pub struct LiftOverResult<'a> {
    pub chromosome: &'a Chromosome,
    pub position: u64,
    pub chain_index: usize,
    pub strand: Strand,
}

#[derive(Debug)]
pub struct PositionLiftOver {
    original_interval: HashMap<String, IntervalTree<u64, TargetRegion>>,
    chain_file: ChainFile,
}

impl PositionLiftOver {
    pub fn new(chain_file: ChainFile) -> PositionLiftOver {
        let mut original_interval = HashMap::new();

        for (i, one_chain) in chain_file.chain_list.iter().enumerate() {
            trace!(
                "Register chain: {}  ({} -> {})",
                one_chain.chain_id,
                one_chain.original_chromosome.name,
                one_chain.new_chromosome.name
            );

            let one_original_interval = original_interval
                .entry(one_chain.original_chromosome.name.clone())
                .or_insert_with(IntervalTree::new);
            // let one_query_interval = query_interval
            //     .entry(one_chain.query_chromosome.name.clone())
            //     .or_insert_with(IntervalTree::new);

            chain_register_interval_tree(
                one_chain,
                one_original_interval,
                // one_query_interval,
                chain_file.original_chromosome_name_to_index[&one_chain.original_chromosome.name],
                chain_file.new_chromosome_name_to_index[&one_chain.new_chromosome.name],
                i,
            );
        }

        PositionLiftOver {
            original_interval,
            //query_interval,
            chain_file,
        }
    }

    pub fn load<R: Read>(chain_reader: R) -> Result<PositionLiftOver, LiftOverError> {
        let chain_file = ChainFile::load(chain_reader)?;
        Ok(PositionLiftOver::new(chain_file))
    }

    fn get_original_interval(&self, chromosome: &str) -> Option<&IntervalTree<u64, TargetRegion>> {
        if let Some(val) = self.original_interval.get(chromosome) {
            Some(val)
        } else if let Some(val) = self.original_interval.get(&format!("chr{}", chromosome)) {
            Some(val)
        } else if chromosome.len() >= 4 && chromosome.starts_with("chr") {
            if let Some(val) = self.original_interval.get(&chromosome[3..]) {
                Some(val)
            } else {
                None
            }
        } else {
            None
        }
    }
    pub fn chain_list(&self) -> &[Chain] {
        &self.chain_file.chain_list
    }

    pub fn original_chromosomes(&self) -> &[Chromosome] {
        &self.chain_file.original_chromosomes
    }

    pub fn new_chromosomes(&self) -> &[Chromosome] {
        &self.chain_file.new_chromosomes
    }

    pub fn original_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.chain_file.original_chromosome_by_name(name)
    }

    pub fn new_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.chain_file.new_chromosome_by_name(name)
    }

    // LiftOver a genomic coordinate. A position should be zero-based.
    pub fn lift_position(&self, chromosome: &str, position: u64) -> Vec<LiftOverResult> {
        let mut results = Vec::new();
        for loop_one in self.search_target(chromosome, position..(position + 1)) {
            for data in loop_one {
                // skip indel
                if (data.original_end - data.original_start) != (data.new_end - data.new_start) {
                    continue;
                }

                if data.is_in_gap {
                    continue;
                }

                results.push(LiftOverResult {
                    chromosome: self
                        .new_chromosomes()
                        .get(data.new_chromosome_index)
                        .unwrap(),
                    position: match data.strand {
                        Strand::Forward => data.new_start + (position - data.original_start),
                        Strand::Reverse => data.new_start + (data.original_end - position) - 1,
                    },
                    chain_index: data.chain_index,
                    strand: data.strand,
                });
            }
        }
        results
    }

    pub fn search_target<'a>(
        &'a self,
        chromosome: &str,
        range: Range<u64>,
    ) -> Vec<Vec<&'a TargetRegion>> {
        let mut groups = HashMap::new();
        if let Some(interval_tree) = self.get_original_interval(chromosome) {
            for one in interval_tree.find(range.clone()) {
                let data: &TargetRegion = one.data();
                if data.original_end == data.original_start && range.start == data.original_start {
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
                Strand::Forward => (x.original_start, x.original_end, x.new_start, x.new_end),
                Strand::Reverse => (
                    x.original_start,
                    x.original_end,
                    u64::MAX - x.new_start,
                    u64::MAX - x.new_end,
                ),
            });
        }

        groups.into_iter().map(|x| x.1).collect()
    }

    pub fn lift_region(&self, chromosome: &str, range: Range<u64>) -> Vec<LiftRegionResult> {
        let mut result = Vec::new();
        for one_chain in self.search_target(chromosome, range.clone()) {
            // println!("lift region: {:?}", one_chain);
            if one_chain.len() == 1 {
                let one_target = one_chain[0];
                // println!(
                //     "check one target {} {} {}",
                //     one_target.contains_reference_position(start),
                //     one_target.contains_reference_position(end),
                //     !one_target.is_indel()
                // );
                if one_target.original_start <= range.start
                    && range.end <= one_target.original_end
                    && !one_target.is_in_gap
                {
                    let start_offset = range.start - one_target.original_start;
                    let end_offset = one_target.original_end - range.end;

                    result.push(LiftRegionResult {
                        chromosome: &self.new_chromosomes()[one_chain[0].new_chromosome_index],
                        start: match one_target.strand {
                            Strand::Forward => one_target.new_start + start_offset,
                            Strand::Reverse => one_target.new_start + end_offset,
                        },
                        end: match one_target.strand {
                            Strand::Forward => one_target.new_end - end_offset,
                            Strand::Reverse => one_target.new_end - start_offset,
                        },
                        strand: one_target.strand,
                        changes: vec![RegionChangeOp::Aligned(range.end - range.start)],
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
                if one_chain[0].original_start > range.start {
                    //println!("skipped first target");
                    continue; // out of chain
                }

                let first_offset = range.start - one_chain[0].original_start;
                let query_start = if !one_chain[0].is_in_gap {
                    // aligned region
                    changes.push(RegionChangeOp::Aligned(
                        one_chain[0].original_end - one_chain[0].original_start - first_offset,
                    ));

                    match one_chain[0].strand {
                        Strand::Forward => one_chain[0].new_start + first_offset,
                        Strand::Reverse => one_chain[0].new_end - first_offset,
                    }
                } else {
                    // INDEL
                    changes.push(RegionChangeOp::Deletion(
                        one_chain[0].original_end - one_chain[0].original_start - first_offset,
                    ));

                    if one_chain[1].is_in_gap {
                        // INDEL should not be repeated
                        unreachable!()
                    }

                    match one_chain[0].strand {
                        Strand::Forward => one_chain[0].new_end,
                        Strand::Reverse => one_chain[0].new_start,
                    }
                };

                for one_target in &one_chain[1..(one_chain.len() - 1)] {
                    if !one_target.is_in_gap {
                        // aligned region
                        changes.push(RegionChangeOp::Aligned(
                            one_target.original_end - one_target.original_start,
                        ));
                    } else {
                        let insert_length = one_target.new_end - one_target.new_start;
                        let delete_length = one_target.original_end - one_target.original_start;
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
                if last_target.original_end < range.end {
                    // out of chain
                    //println!("skipped last target");
                    continue;
                }
                let last_offset = last_target.original_end - range.end;
                let query_end = if !last_target.is_indel() {
                    // aligned region
                    changes.push(RegionChangeOp::Aligned(
                        last_target.original_end - last_target.original_start - last_offset,
                    ));

                    match one_chain[0].strand {
                        Strand::Forward => last_target.new_end - last_offset,
                        Strand::Reverse => last_target.new_start + last_offset,
                    }
                } else {
                    // INDEL
                    changes.push(RegionChangeOp::Deletion(
                        last_target.original_end - last_target.original_start - last_offset,
                    ));

                    match one_chain[0].strand {
                        Strand::Forward => last_target.new_start,
                        Strand::Reverse => last_target.new_end,
                    }
                };

                result.push(LiftRegionResult {
                    chromosome: &self.new_chromosomes()[one_chain[0].new_chromosome_index],
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
    original_interval: &mut IntervalTree<u64, TargetRegion>,
    //query_interval: &mut IntervalTree<u64, TargetRegion>,
    original_chromosome_index: usize,
    new_chromosome_index: usize,
    chain_index: usize,
) {
    assert_eq!(chain.original_strand, Strand::Forward);

    let mut original_current = chain.original_start;
    let mut new_current = chain.new_start;
    for one in chain.chain_interval.iter() {
        // insert un-gapped region
        let original_next = original_current + one.size;
        let new_next = new_current + one.size;

        register_one_interval(
            chain,
            original_interval,
            //query_interval,
            original_chromosome_index,
            new_chromosome_index,
            chain_index,
            original_current,
            original_next,
            new_current,
            new_next,
            false,
        );

        original_current = original_next;
        new_current = new_next;

        // insert gapped region
        if let Some(diff_ref) = one.difference_original {
            if let Some(diff_query) = one.difference_new {
                let original_next = original_current + diff_ref;
                let new_next = new_current + diff_query;

                register_one_interval(
                    chain,
                    original_interval,
                    //query_interval,
                    original_chromosome_index,
                    new_chromosome_index,
                    chain_index,
                    original_current,
                    original_next,
                    new_current,
                    new_next,
                    true,
                );

                original_current = original_next;
                new_current = new_next;
            } else {
                panic!("Invalid format")
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn register_one_interval(
    chain: &Chain,
    original_interval: &mut IntervalTree<u64, TargetRegion>,
    //query_interval: &mut IntervalTree<u64, TargetRegion>,
    original_chromosome_index: usize,
    new_chromosome_index: usize,
    chain_index: usize,
    original_current: u64,
    original_next: u64,
    new_current: u64,
    new_next: u64,
    is_in_gap: bool,
) {
    match chain.new_strand {
        Strand::Forward => {
            let target_region = TargetRegion {
                original_chromosome_index,
                new_chromosome_index,
                chain_index,
                original_start: original_current,
                original_end: original_next,
                new_start: new_current,
                new_end: new_next,
                strand: Strand::Forward,
                is_in_gap,
            };

            let register_original_next = if original_current == original_next {
                original_next + 1
            } else {
                original_next
            };
            // let register_query_next = if query_current == query_next {
            //     query_next + 1
            // } else {
            //     query_next
            // };

            original_interval.insert(
                original_current..register_original_next,
                target_region.clone(),
            );
            //query_interval.insert(query_current..register_query_next, target_region);
        }
        Strand::Reverse => {
            let target_region = TargetRegion {
                original_chromosome_index,
                new_chromosome_index,
                chain_index,
                original_start: original_current,
                original_end: original_next,
                new_start: chain.new_chromosome.length - new_next,
                new_end: chain.new_chromosome.length - new_current,
                strand: Strand::Reverse,
                is_in_gap,
            };

            let register_original_next = if original_current == original_next {
                original_next + 1
            } else {
                original_next
            };
            // let register_query_next = if query_current == query_next {
            //     query_next + 1
            // } else {
            //     query_next
            // };

            original_interval.insert(
                original_current..register_original_next,
                target_region.clone(),
            );
            // query_interval.insert(
            //     (chain.query_chromosome.length - register_query_next)
            //         ..(chain.query_chromosome.length - query_current),
            //     target_region,
            // );
        }
    }
}

#[cfg(test)]
mod test;

#[cfg(all(test, feature = "test_all"))]
mod test2;
