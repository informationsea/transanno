use crate::{reverse_complement, GenomeSequence, LiftOverError, Variant};
use regex::Regex;

use log::{info, trace};
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::io::{BufRead, BufReader, Read, Write};
use std::str::FromStr;
use std::u64;

lazy_static! {
    static ref SPACE: Regex = Regex::new("\\s").unwrap();
}

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
    pub difference_reference: Option<u64>,
    pub difference_query: Option<u64>,
}

impl ChainInterval {
    pub fn is_valid(&self) -> bool {
        self.size != 0
            && (self.difference_query.unwrap_or(0) != 0
                || self.difference_reference.unwrap_or(0) != 0)
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Chain {
    pub chain_interval: Vec<ChainInterval>,
    pub score: i64,
    pub reference_chromosome: Chromosome,
    pub reference_strand: Strand,
    pub reference_start: u64,
    pub reference_end: u64,
    pub query_strand: Strand,
    pub query_chromosome: Chromosome,
    pub query_start: u64,
    pub query_end: u64,
    pub chain_id: String,
}

impl Chain {
    fn cleanup(chain_interval: &[ChainInterval]) -> Vec<ChainInterval> {
        if chain_interval.is_empty() {
            return Vec::new();
        }

        chain_interval.iter().fold(Vec::new(), |mut acc, x| {
            if let Some(last) = acc.pop() {
                if last.difference_query.unwrap_or(0) == 0
                    && last.difference_reference.unwrap_or(0) == 0
                {
                    acc.push(ChainInterval {
                        size: last.size + x.size,
                        difference_query: x.difference_query,
                        difference_reference: x.difference_reference,
                    });
                } else if x.size == 0 {
                    acc.push(ChainInterval {
                        size: last.size,
                        difference_query: Some(
                            last.difference_query.unwrap_or(0) + x.difference_query.unwrap_or(0),
                        ),
                        difference_reference: Some(
                            last.difference_reference.unwrap_or(0)
                                + x.difference_reference.unwrap_or(0),
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

    pub fn left_align<G: GenomeSequence>(
        &self,
        reference: &mut G,
        query: &mut G,
    ) -> Result<Chain, LiftOverError> {
        if self.reference_strand == Strand::Reverse {
            return Err(LiftOverError::ReferenceStrandShouldForward);
        }

        let mut new_intervals = Vec::new();

        let mut current_reference = self.reference_start;
        let mut current_query = self.query_start;
        let mut remain_size = 0;

        for one_interval in self.chain_interval.iter() {
            if one_interval.difference_query.unwrap_or(0) == 1
                && one_interval.difference_reference.unwrap_or(0) == 1
            {
                // println!("concatenate {} {}", one_interval.size, remain_size);
                remain_size += one_interval.size + 1;
                continue;
            }

            let next_reference = current_reference + one_interval.size + remain_size;
            let next_query = current_query + one_interval.size + remain_size;
            trace!(
                "next reference: {} - {} {} {}",
                next_reference,
                current_reference,
                one_interval.size,
                remain_size
            );

            let reference_seq = reference.get_sequence(
                &self.reference_chromosome.name,
                next_reference,
                next_reference + one_interval.difference_reference.unwrap_or(0),
            )?;
            let query_seq = match self.query_strand {
                Strand::Forward => query.get_sequence(
                    &self.query_chromosome.name,
                    next_query,
                    next_query + one_interval.difference_query.unwrap_or(0),
                )?,
                Strand::Reverse => reverse_complement(&query.get_sequence(
                    &self.query_chromosome.name,
                    self.query_chromosome.length
                        - (next_query + one_interval.difference_query.unwrap_or(0)),
                    self.query_chromosome.length - next_query,
                )?),
            };
            let variant = Variant {
                chromosome: self.reference_chromosome.name.to_string(),
                position: next_reference,
                reference: reference_seq,
                alternative: vec![query_seq],
            };
            //println!("before normalization variant: {:?}", variant);
            let normalized = variant
                .normalize(reference)?
                .truncate_left_most_nucleotide_if_allele_starts_with_same();
            //println!(" after normalization variant: {:?}", normalized);

            let offset = if next_reference > normalized.position {
                next_reference - normalized.position
            } else {
                0
            };

            let mut offset_updated = if offset == 0 {
                0
            } else {
                let offset_reference_seq = reference.get_sequence(
                    &self.reference_chromosome.name,
                    next_reference - offset,
                    next_reference,
                )?;
                let offset_query_seq = match self.query_strand {
                    Strand::Forward => query.get_sequence(
                        &self.query_chromosome.name,
                        next_query - offset,
                        next_query,
                    )?,
                    Strand::Reverse => reverse_complement(&query.get_sequence(
                        &self.query_chromosome.name,
                        self.query_chromosome.length - next_query,
                        self.query_chromosome.length - (next_query - offset),
                    )?),
                };
                if offset_reference_seq == offset_query_seq {
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
                difference_reference: one_interval
                    .difference_reference
                    .map(|_| normalized.reference.len() as u64),
                difference_query: one_interval
                    .difference_query
                    .map(|_| normalized.alternative[0].len() as u64),
            };

            let original_current_reference = current_reference
                + one_interval.size
                + one_interval.difference_reference.unwrap_or(0)
                + remain_size;
            current_reference +=
                one_new_interval.size + one_new_interval.difference_reference.unwrap_or(0);
            current_query += one_new_interval.size + one_new_interval.difference_query.unwrap_or(0);
            // println!(
            //     "original: {}  / current: {}",
            //     original_current_reference, current_reference
            // );
            remain_size = original_current_reference - current_reference;
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
    pub reference_chromosomes: Vec<Chromosome>,
    pub query_chromosomes: Vec<Chromosome>,

    pub reference_chromosome_name_to_index: HashMap<String, usize>,
    pub query_chromosome_name_to_index: HashMap<String, usize>,
}

impl ChainFile {
    pub fn left_align<G: GenomeSequence>(
        &self,
        reference: &mut G,
        query: &mut G,
    ) -> Result<ChainFile, LiftOverError> {
        info!("start  chain file left align");
        let mut new_chain_list = Vec::new();
        for one_chain in self.chain_list.iter() {
            new_chain_list.push(one_chain.left_align(reference, query)?);
        }
        info!("finish chain file left align");

        Ok(ChainFile {
            chain_list: new_chain_list,
            reference_chromosomes: self.reference_chromosomes.clone(),
            query_chromosomes: self.query_chromosomes.clone(),
            reference_chromosome_name_to_index: self.reference_chromosome_name_to_index.clone(),
            query_chromosome_name_to_index: self.query_chromosome_name_to_index.clone(),
        })
    }

    pub fn load<R: Read>(chain_file: R) -> Result<ChainFile, LiftOverError> {
        let mut reader = BufReader::new(chain_file);
        let mut status = LiftOverReadStatus::Outside;
        let mut line = String::new();
        let mut line_num = 0;

        let mut chain_list = Vec::new();
        let mut current_chain: Option<Chain> = None;

        let mut reference_chromosomes: Vec<Chromosome> = Vec::new();
        let mut query_chromosomes: Vec<Chromosome> = Vec::new();

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
                        reference_chromosome: Chromosome {
                            name: elements[2].to_string(),
                            length: elements[3].parse()?,
                        },
                        reference_strand: elements[4].parse()?,
                        reference_start: elements[5].parse()?,
                        reference_end: elements[6].parse()?,
                        query_chromosome: Chromosome {
                            name: elements[7].to_string(),
                            length: elements[8].parse()?,
                        },
                        query_strand: elements[9].parse()?,
                        query_start: elements[10].parse()?,
                        query_end: elements[11].parse()?,
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
                                difference_reference: Some(elements[1].parse()?),
                                difference_query: Some(elements[2].parse()?),
                            });
                    } else if elements.len() == 1 {
                        if let Some(mut current_chain) = current_chain.take() {
                            current_chain.chain_interval.push(ChainInterval {
                                size: elements[0].parse()?,
                                difference_reference: None,
                                difference_query: None,
                            });

                            // last entry
                            status = LiftOverReadStatus::Outside;

                            // register chromosome
                            if let Some(registered_chromosome) = reference_chromosomes
                                .iter()
                                .find(|x| x.name == current_chain.reference_chromosome.name)
                            {
                                if registered_chromosome.length
                                    != current_chain.reference_chromosome.length
                                {
                                    return Err(LiftOverError::InvalidChromosomeLength(line_num));
                                }
                            } else {
                                reference_chromosomes
                                    .push(current_chain.reference_chromosome.clone());
                            }

                            if let Some(registered_chromosome) = query_chromosomes
                                .iter()
                                .find(|x| x.name == current_chain.query_chromosome.name)
                            {
                                if registered_chromosome.length
                                    != current_chain.query_chromosome.length
                                {
                                    return Err(LiftOverError::InvalidChromosomeLength(line_num));
                                }
                            } else {
                                query_chromosomes.push(current_chain.query_chromosome.clone());
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
            reference_chromosome_name_to_index: reference_chromosomes
                .iter()
                .enumerate()
                .map(|(i, x)| (x.name.clone(), i))
                .collect(),
            query_chromosome_name_to_index: query_chromosomes
                .iter()
                .enumerate()
                .map(|(i, x)| (x.name.clone(), i))
                .collect(),
            reference_chromosomes,
            query_chromosomes,
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
                one_chain.reference_chromosome.name,
                one_chain.reference_chromosome.length,
                one_chain.reference_strand,
                one_chain.reference_start,
                one_chain.reference_end,
                one_chain.query_chromosome.name,
                one_chain.query_chromosome.length,
                one_chain.query_strand,
                one_chain.query_start,
                one_chain.query_end,
                one_chain.chain_id
            )?;
            for one_interval in one_chain.chain_interval.iter() {
                if let (Some(diff_ref), Some(diff_query)) = (
                    one_interval.difference_reference,
                    one_interval.difference_query,
                ) {
                    writeln!(
                        writer,
                        "{}\t{}\t{}",
                        one_interval.size, diff_ref, diff_query
                    )?;
                } else {
                    writeln!(writer, "{}", one_interval.size)?;
                }
            }
        }
        Ok(())
    }

    pub fn query_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.query_chromosome_name_to_index
            .get(name)
            .map(|x| &self.query_chromosomes[*x])
    }

    pub fn reference_chromosome_by_name(&self, name: &str) -> Option<&Chromosome> {
        self.reference_chromosome_name_to_index
            .get(name)
            .map(|x| &self.reference_chromosomes[*x])
    }
}

#[cfg(test)]
mod test;
