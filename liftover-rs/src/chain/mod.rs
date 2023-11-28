use crate::{reverse_complement, GenomeSequence, LiftOverError, Variant};
use once_cell::sync::Lazy;
use regex::Regex;

use log::{error, info, trace, warn};
use std::collections::{HashMap, HashSet};
use std::fmt::{self, Display};
use std::io::{BufRead, BufReader, Read, Write};
use std::str::FromStr;
use std::u64;

static SPACE: Lazy<Regex> = Lazy::new(|| Regex::new("\\s").unwrap());

#[derive(Debug, Hash, Default, Clone, PartialEq, Eq)]
pub struct Chromosome {
    pub name: String,
    pub length: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy, PartialOrd, Ord)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Default for Strand {
    fn default() -> Self {
        Strand::Forward
    }
}

impl FromStr for Strand {
    type Err = LiftOverError;
    fn from_str(s: &str) -> Result<Self, LiftOverError> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => Err(LiftOverError::ParseStrandError),
        }
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}

#[derive(Debug, PartialEq)]
enum LiftOverReadStatus {
    Outside,
    InChain,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChainInterval {
    pub size: u64,
    pub difference_original: Option<u64>,
    pub difference_new: Option<u64>,
}

impl ChainInterval {
    pub fn is_valid(&self) -> bool {
        self.size != 0
            && (self.difference_new.unwrap_or(0) != 0 || self.difference_original.unwrap_or(0) != 0)
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Chain {
    pub chain_interval: Vec<ChainInterval>,
    pub score: i64,
    pub original_chromosome: Chromosome,
    pub original_strand: Strand,
    pub original_start: u64,
    pub original_end: u64,
    pub new_strand: Strand,
    pub new_chromosome: Chromosome,
    pub new_start: u64,
    pub new_end: u64,
    pub chain_id: String,
}

impl Chain {
    fn cleanup(chain_interval: &[ChainInterval]) -> Vec<ChainInterval> {
        if chain_interval.is_empty() {
            return Vec::new();
        }

        chain_interval.iter().fold(Vec::new(), |mut acc, x| {
            if let Some(last) = acc.pop() {
                if last.difference_new.unwrap_or(0) == 0
                    && last.difference_original.unwrap_or(0) == 0
                {
                    acc.push(ChainInterval {
                        size: last.size + x.size,
                        difference_new: x.difference_new,
                        difference_original: x.difference_original,
                    });
                } else if x.size == 0 {
                    acc.push(ChainInterval {
                        size: last.size,
                        difference_new: Some(
                            last.difference_new.unwrap_or(0) + x.difference_new.unwrap_or(0),
                        ),
                        difference_original: Some(
                            last.difference_original.unwrap_or(0)
                                + x.difference_original.unwrap_or(0),
                        ),
                    });
                } else {
                    acc.push(last);
                    acc.push(x.clone());
                }
            } else {
                acc.push(x.clone());
            }
            acc
        })
    }

    pub fn check_sequence_consistency<G1: GenomeSequence, G2: GenomeSequence>(
        &self,
        original_sequence: &mut G1,
        new_sequence: &mut G2,
    ) -> Result<(), LiftOverError> {
        // Check chromosome length
        if let Some(expected_len) = original_sequence
            .get_contig_list()
            .iter()
            .filter(|x| x.0 == self.original_chromosome.name)
            .map(|x| x.1)
            .next()
        {
            if expected_len != self.original_chromosome.length {
                error!(
                    "Length of {} in chain file is {}, but {} in original FASTA is {}",
                    self.original_chromosome.name,
                    self.original_chromosome.length,
                    self.original_chromosome.name,
                    expected_len
                );
                return Err(LiftOverError::UnmatchedOriginalChromosomeLength(
                    self.original_chromosome.name.to_string(),
                    self.original_chromosome.length,
                    expected_len,
                ));
            }
        } else {
            warn!(
                "Chromosome {} is not found in original FASTA",
                self.original_chromosome.name
            );
            return Err(LiftOverError::ChromosomeNotFound(
                self.original_chromosome.name.to_string(),
            ));
        }

        if let Some(expected_len) = new_sequence
            .get_contig_list()
            .iter()
            .filter(|x| x.0 == self.new_chromosome.name)
            .map(|x| x.1)
            .next()
        {
            if expected_len != self.new_chromosome.length {
                error!(
                    "Length of {} in chain file is {}, but {} in new FASTA is {}",
                    self.new_chromosome.name,
                    self.new_chromosome.length,
                    self.new_chromosome.name,
                    expected_len
                );

                return Err(LiftOverError::UnmatchedNewChromosomeLength(
                    self.original_chromosome.name.to_string(),
                    self.original_chromosome.length,
                    expected_len,
                ));
            }
        } else {
            warn!(
                "Chromosome {} is not found in new FASTA",
                self.new_chromosome.name
            );
            return Err(LiftOverError::ChromosomeNotFound(
                self.new_chromosome.name.to_string(),
            ));
        }
        Ok(())
    }

    pub fn left_align<G: GenomeSequence>(
        &self,
        original_sequence: &mut G,
        new_sequence: &mut G,
    ) -> Result<Chain, LiftOverError> {
        if self.original_strand == Strand::Reverse {
            return Err(LiftOverError::OriginalStrandShouldForward);
        }
        self.check_sequence_consistency(original_sequence, new_sequence)?;

        let mut new_intervals = Vec::new();

        let mut current_original = self.original_start;
        let mut current_new = self.new_start;
        let mut remain_size = 0;

        for (i, one_interval) in self.chain_interval.iter().enumerate() {
            trace!(
                "working {}/{} - {:?} (chain-id:{})",
                i,
                self.chain_interval.len(),
                one_interval,
                self.chain_id
            );
            if one_interval.difference_new.unwrap_or(0) == 1
                && one_interval.difference_original.unwrap_or(0) == 1
            {
                // println!("concatenate {} {}", one_interval.size, remain_size);
                remain_size += one_interval.size + 1;
                continue;
            }

            let next_original = current_original + one_interval.size + remain_size;
            let next_new = current_new + one_interval.size + remain_size;
            trace!(
                "next reference: {}:{} - {} {} {}",
                &self.original_chromosome.name,
                next_original,
                current_original,
                one_interval.size,
                remain_size
            );

            let original_seq = original_sequence.get_sequence(
                &self.original_chromosome.name,
                next_original,
                next_original + one_interval.difference_original.unwrap_or(0),
            )?;

            match self.new_strand {
                Strand::Forward => {
                    trace!(
                        "next new+: {}:{} - {} {} {}",
                        &self.new_chromosome.name,
                        next_new,
                        next_new + one_interval.difference_new.unwrap_or(0),
                        one_interval.size,
                        remain_size
                    );
                }
                Strand::Reverse => {
                    trace!(
                        "next new-: {}:{} - {} {} {}",
                        &self.new_chromosome.name,
                        self.new_chromosome.length
                            - (next_new + one_interval.difference_new.unwrap_or(0)),
                        self.new_chromosome.length - next_new,
                        one_interval.size,
                        remain_size
                    );
                }
            };

            let new_seq = match self.new_strand {
                Strand::Forward => new_sequence.get_sequence(
                    &self.new_chromosome.name,
                    next_new,
                    next_new + one_interval.difference_new.unwrap_or(0),
                )?,
                Strand::Reverse => reverse_complement(&new_sequence.get_sequence(
                    &self.new_chromosome.name,
                    self.new_chromosome.length
                        - (next_new + one_interval.difference_new.unwrap_or(0)),
                    self.new_chromosome.length - next_new,
                )?),
            };
            trace!("new seq ok");

            let do_not_normalize = new_seq.contains(&b'N') || original_seq.contains(&b'N');

            let variant = Variant {
                chromosome: self.original_chromosome.name.to_string(),
                position: next_original,
                reference: original_seq,
                alternative: vec![new_seq],
            };
            //println!("before normalization variant: {:?}", variant);

            let normalized = if do_not_normalize {
                variant
            } else {
                variant
                    .normalize(original_sequence)?
                    .truncate_left_most_nucleotide_if_allele_starts_with_same()
            };
            //println!(" after normalization variant: {:?}", normalized);

            let offset = if next_original > normalized.position {
                next_original - normalized.position
            } else {
                0
            };

            let mut offset_updated = if offset == 0 {
                0
            } else {
                trace!("offset fetch: {}/{}/{}", next_original, next_new, offset);
                let offset_original_seq = original_sequence.get_sequence(
                    &self.original_chromosome.name,
                    next_original - offset,
                    next_original,
                )?;
                let offset_query_seq = match self.new_strand {
                    Strand::Forward => new_sequence.get_sequence(
                        &self.new_chromosome.name,
                        next_new - offset,
                        next_new,
                    )?,
                    Strand::Reverse => reverse_complement(&new_sequence.get_sequence(
                        &self.new_chromosome.name,
                        self.new_chromosome.length - next_new,
                        self.new_chromosome.length - (next_new - offset),
                    )?),
                };
                if offset_original_seq == offset_query_seq {
                    offset
                } else {
                    0
                }
            };

            if offset_updated > one_interval.size + remain_size {
                offset_updated = one_interval.size + remain_size;
            }

            trace!(
                "interval size: {} {} {}",
                one_interval.size,
                offset_updated,
                remain_size
            );

            let one_new_interval = ChainInterval {
                size: one_interval.size + remain_size - offset_updated,
                difference_original: one_interval
                    .difference_original
                    .map(|_| normalized.reference.len() as u64),
                difference_new: one_interval
                    .difference_new
                    .map(|_| normalized.alternative[0].len() as u64),
            };

            let original_current_original = current_original
                + one_interval.size
                + one_interval.difference_original.unwrap_or(0)
                + remain_size;
            current_original +=
                one_new_interval.size + one_new_interval.difference_original.unwrap_or(0);
            current_new += one_new_interval.size + one_new_interval.difference_new.unwrap_or(0);
            // println!(
            //     "original: {}  / current: {}",
            //     original_current_reference, current_reference
            // );
            remain_size = original_current_original - current_original;
            new_intervals.push(one_new_interval);
        }

        Ok(Chain {
            chain_interval: Chain::cleanup(&new_intervals),
            ..self.clone()
        })
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct ChainFile {
    pub chain_list: Vec<Chain>,
    pub original_chromosomes: Vec<Chromosome>,
    pub new_chromosomes: Vec<Chromosome>,

    pub original_chromosome_name_to_index: HashMap<String, usize>,
    pub new_chromosome_name_to_index: HashMap<String, usize>,
}

impl ChainFile {
    pub fn left_align<G: GenomeSequence>(
        &self,
        original_sequence: &mut G,
        new_sequence: &mut G,
    ) -> Result<ChainFile, LiftOverError> {
        info!("start chain file left align");
        let mut new_chain_list = Vec::new();
        let mut warned_original_chrom = HashSet::new();
        let mut warned_new_chrom = HashSet::new();
        for one_chain in self.chain_list.iter() {
            if let Some(length) =
                original_sequence.get_contig_length(&one_chain.original_chromosome.name)
            {
                if one_chain.original_chromosome.length != length {
                    if !warned_original_chrom.contains(&one_chain.original_chromosome.name) {
                        warn!("Length of chromosome {} is not match between chain file and original FASTA. Skipping...", one_chain.original_chromosome.name);
                        warned_original_chrom.insert(one_chain.original_chromosome.name.clone());
                    }
                    continue;
                }
            } else {
                if !warned_original_chrom.contains(&one_chain.original_chromosome.name) {
                    warn!(
                        "Chromosome {} is not found in original FASTA. Skipping...",
                        one_chain.original_chromosome.name
                    );
                    warned_original_chrom.insert(one_chain.original_chromosome.name.clone());
                }
                continue;
            }

            if let Some(length) = new_sequence.get_contig_length(&one_chain.new_chromosome.name) {
                if one_chain.new_chromosome.length != length {
                    if !warned_new_chrom.contains(&one_chain.new_chromosome.name) {
                        warn!("Length of chromosome {} is not match between chain file and new FASTA. Skipping...", one_chain.new_chromosome.name);
                        warned_new_chrom.insert(one_chain.new_chromosome.name.clone());
                    }
                    continue;
                }
            } else {
                if !warned_new_chrom.contains(&one_chain.new_chromosome.name) {
                    warn!(
                        "Chromosome {} is not found in new FASTA. Skipping...",
                        one_chain.new_chromosome.name
                    );
                    warned_new_chrom.insert(one_chain.new_chromosome.name.clone());
                }
                continue;
            }

            new_chain_list.push(one_chain.left_align(original_sequence, new_sequence)?);
        }
        if new_chain_list.is_empty() {
            error!("No chain found. You may use wrong FASTA.");
            return Err(LiftOverError::NoChainFoundError);
        }
        info!("finish chain file left align");

        Ok(ChainFile {
            chain_list: new_chain_list,
            original_chromosomes: self.original_chromosomes.clone(),
            new_chromosomes: self.new_chromosomes.clone(),
            original_chromosome_name_to_index: self.original_chromosome_name_to_index.clone(),
            new_chromosome_name_to_index: self.new_chromosome_name_to_index.clone(),
        })
    }

    pub fn load<R: Read>(chain_file: R) -> Result<ChainFile, LiftOverError> {
        let mut reader = BufReader::new(chain_file);
        let mut status = LiftOverReadStatus::Outside;
        let mut line = String::new();
        let mut line_num = 0;

        let mut chain_list = Vec::new();
        let mut current_chain: Option<Chain> = None;

        let mut original_chromosomes: Vec<Chromosome> = Vec::new();
        let mut new_chromosomes: Vec<Chromosome> = Vec::new();

        loop {
            line_num += 1;
            line.clear();
            let read_bytes = reader.read_line(&mut line)?;
            if read_bytes == 0 {
                // last line
                break;
            }
            let trim_line = line.trim();

            if trim_line.starts_with('#') {
                // Skip comment line
                continue;
            }
            //println!("line {}({:?}): {}", line_num, status, trim_line);
            match status {
                LiftOverReadStatus::Outside => {
                    if trim_line.is_empty() {
                        continue;
                    }

                    let elements: Vec<&str> = SPACE.split(trim_line).collect();
                    if elements.len() != 13 {
                        return Err(LiftOverError::InvalidNumberOfHeader(line_num));
                    }

                    if elements[0] != "chain" {
                        return Err(LiftOverError::NoChainHeaderFound(line_num));
                    }

                    if elements[4] != "+" {
                        return Err(LiftOverError::InvalidStrand(line_num));
                    }

                    status = LiftOverReadStatus::InChain;
                    current_chain = Some(Chain {
                        score: elements[1].parse()?,
                        original_chromosome: Chromosome {
                            name: elements[2].to_string(),
                            length: elements[3].parse()?,
                        },
                        original_strand: elements[4].parse()?,
                        original_start: elements[5].parse()?,
                        original_end: elements[6].parse()?,
                        new_chromosome: Chromosome {
                            name: elements[7].to_string(),
                            length: elements[8].parse()?,
                        },
                        new_strand: elements[9].parse()?,
                        new_start: elements[10].parse()?,
                        new_end: elements[11].parse()?,
                        chain_id: elements[12].to_string(),
                        chain_interval: Vec::new(),
                    });
                }
                LiftOverReadStatus::InChain => {
                    let elements: Vec<&str> = SPACE.split(trim_line).collect();
                    //println!("in chain: {:?}", elements);
                    if elements.len() == 3 {
                        current_chain
                            .as_mut()
                            .unwrap()
                            .chain_interval
                            .push(ChainInterval {
                                size: elements[0].parse()?,
                                difference_original: Some(elements[1].parse()?),
                                difference_new: Some(elements[2].parse()?),
                            });
                    } else if elements.len() == 1 {
                        if let Some(mut current_chain) = current_chain.take() {
                            current_chain.chain_interval.push(ChainInterval {
                                size: elements[0].parse()?,
                                difference_original: None,
                                difference_new: None,
                            });

                            // last entry
                            status = LiftOverReadStatus::Outside;

                            // register chromosome
                            if let Some(registered_chromosome) = original_chromosomes
                                .iter()
                                .find(|x| x.name == current_chain.original_chromosome.name)
                            {
                                if registered_chromosome.length
                                    != current_chain.original_chromosome.length
                                {
                                    return Err(LiftOverError::InvalidChromosomeLength(line_num));
                                }
                            } else {
                                original_chromosomes
                                    .push(current_chain.original_chromosome.clone());
                            }

                            if let Some(registered_chromosome) = new_chromosomes
                                .iter()
                                .find(|x| x.name == current_chain.new_chromosome.name)
                            {
                                if registered_chromosome.length
                                    != current_chain.new_chromosome.length
                                {
                                    return Err(LiftOverError::InvalidChromosomeLength(line_num));
                                }
                            } else {
                                new_chromosomes.push(current_chain.new_chromosome.clone());
                            }

                            chain_list.push(current_chain);
                        }
                    } else {
                        return Err(LiftOverError::InvalidNumberOfColumns(line_num));
                    }
                }
            }
        }

        Ok(ChainFile {
            chain_list,
            original_chromosome_name_to_index: original_chromosomes
                .iter()
                .enumerate()
                .map(|(i, x)| (x.name.clone(), i))
                .collect(),
            new_chromosome_name_to_index: new_chromosomes
                .iter()
                .enumerate()
                .map(|(i, x)| (x.name.clone(), i))
                .collect(),
            original_chromosomes,
            new_chromosomes,
        })
    }

    pub fn write<W: Write>(&self, writer: &mut W) -> Result<(), LiftOverError> {
        for (i, one_chain) in self.chain_list.iter().enumerate() {
            if i != 0 {
                writer.write_all(b"\n")?;
            }
            writeln!(
                writer,
                "chain {} {} {} {} {} {} {} {} {} {} {} {}",
                one_chain.score,
                one_chain.original_chromosome.name,
                one_chain.original_chromosome.length,
                one_chain.original_strand,
                one_chain.original_start,
                one_chain.original_end,
                one_chain.new_chromosome.name,
                one_chain.new_chromosome.length,
                one_chain.new_strand,
                one_chain.new_start,
                one_chain.new_end,
                one_chain.chain_id
            )?;
            for one_interval in one_chain.chain_interval.iter() {
                if let (Some(diff_original), Some(diff_new)) = (
                    one_interval.difference_original,
                    one_interval.difference_new,
                ) {
                    writeln!(
                        writer,
                        "{}\t{}\t{}",
                        one_interval.size, diff_original, diff_new
                    )?;
                } else {
                    writeln!(writer, "{}", one_interval.size)?;
                }
            }
        }
        Ok(())
    }

    pub fn new_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.new_chromosome_name_to_index
            .get(name)
            .map(|x| &self.new_chromosomes[*x])
    }

    pub fn original_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.original_chromosome_name_to_index
            .get(name)
            .map(|x| &self.original_chromosomes[*x])
    }
}

#[cfg(test)]
mod test;

#[cfg(test)]
#[cfg(feature = "test_all")]
mod test2;
