use crate::{GenomeSequence, LiftOverError, LiftOverErrorKind, Variant};

fn extend_to_left_if_empty_allele_exists<G: GenomeSequence>(
    chromosome: &str,
    genome: &mut G,
    position: &mut u64,
    alleles: &mut Vec<Vec<u8>>,
) -> Result<bool, LiftOverError> {
    if alleles.iter().any(|x| x.is_empty()) && *position > 0 {
        *position -= 1;
        let extend_seq = genome.get_sequence(chromosome, *position, *position + 1)?;
        for one in alleles.iter_mut() {
            if one != b"*" {
                one.insert(0, extend_seq[0]);
            }
        }
        Ok(true)
    } else {
        Ok(false)
    }
}

fn truncate_right_most_nucleotide_if_allele_ends_with_same(
    position: u64,
    alleles: &mut Vec<Vec<u8>>,
) -> bool {
    if alleles.len() <= 1
        || alleles.iter().filter(|x| x != &b"*").any(|x| x.is_empty())
        || position == 0
    {
        false
    } else {
        let equal_length = alleles
            .iter()
            .skip(1)
            .filter(|x| x != &b"*")
            .map(|x| {
                x.iter()
                    .rev()
                    .zip(alleles[0].iter().rev())
                    .take_while(|(a, b)| a == b)
                    .count()
            })
            .min();
        if let Some(equal_length) = equal_length {
            if equal_length > 0 {
                for one_allele in alleles.iter_mut() {
                    if one_allele != b"*" {
                        one_allele.truncate(one_allele.len() - equal_length);
                    }
                }
                return true;
            }
        }
        false
    }
}

impl Variant {
    pub fn normalize<G: GenomeSequence>(&self, genome: &mut G) -> Result<Variant, LiftOverError> {
        // do not normalize if all alleles are equal
        if self.alternative.iter().all(|x| x == &self.reference) {
            return Ok(self.clone());
        }
        // skip normalize for SNP
        if self.alternative.iter().all(|x| x.len() == 1) && self.reference.len() == 1 {
            return Ok(self.clone());
        }
        let reference_seq = genome.get_sequence(
            &self.chromosome,
            self.position,
            self.position + self.reference.len() as u64,
        )?;
        if reference_seq != self.reference {
            return Err(LiftOverErrorKind::DifferentReference.into());
        }
        let mut position = self.position;
        let mut alleles = Vec::new();
        alleles.push(reference_seq);
        for one in &self.alternative {
            alleles.push(one.to_vec());
        }

        extend_to_left_if_empty_allele_exists(
            &self.chromosome,
            genome,
            &mut position,
            &mut alleles,
        )?;
        // doi: 10.1093/bioinformatics/btv112
        loop {
            //println!("looping 1 {:?} {:?}", position, alleles);
            let truncated =
                truncate_right_most_nucleotide_if_allele_ends_with_same(position, &mut alleles);
            let extended = extend_to_left_if_empty_allele_exists(
                &self.chromosome,
                genome,
                &mut position,
                &mut alleles,
            )?;
            if !truncated && !extended {
                break;
            }
        }

        let common_prefix_length = (0..(alleles.iter().map(|x| x.len()).min().unwrap() - 1))
            .take_while(|i| alleles.iter().skip(1).all(|x| x[*i] == alleles[0][*i]))
            .count();

        if common_prefix_length > 0 {
            for one_allele in alleles.iter_mut() {
                one_allele.drain(0..(common_prefix_length));
            }
            position += common_prefix_length as u64;
        }

        Ok(Variant {
            chromosome: self.chromosome.clone(),
            position,
            reference: alleles.remove(0),
            alternative: alleles,
        })
    }

    // allow empty alleles
    pub fn truncate_left_most_nucleotide_if_allele_starts_with_same(&self) -> Variant {
        let minimum_length = self
            .alternative
            .iter()
            .filter(|x| x != &b"*")
            .map(|x| x.len())
            .chain(vec![self.reference.len()])
            .min()
            .unwrap();
        let mut matched_length = 0;
        for i in 0..minimum_length {
            if self
                .alternative
                .iter()
                .filter(|x| x != &b"*")
                .all(|x| x[i] == self.reference[i])
            {
                matched_length = i + 1;
            } else {
                break;
            }
        }

        Variant {
            chromosome: self.chromosome.clone(),
            position: self.position + matched_length as u64,
            reference: self.reference[matched_length..].to_vec(),
            alternative: self
                .alternative
                .iter()
                .map(|x| {
                    if x != b"*" {
                        x[matched_length..].to_vec()
                    } else {
                        b"*".to_vec()
                    }
                })
                .collect(),
        }
    }
}

#[cfg(test)]
mod test;
