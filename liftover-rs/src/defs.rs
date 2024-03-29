use super::LiftOverError;
use crate::vcfparse::VCFRecord;
use bio::io::fasta::IndexedReader;
//use log::trace;
use std::fmt;
use std::io::{Read, Seek};
use std::str;

pub trait GenomeSequence: std::fmt::Debug {
    fn sequence(
        &mut self,
        chromosome: &str,
        start: u64,
        stop: u64,
        text: &mut Vec<u8>,
    ) -> Result<(), LiftOverError>;

    fn get_sequence(
        &mut self,
        chromosome: &str,
        start: u64,
        stop: u64,
    ) -> Result<Vec<u8>, LiftOverError> {
        let mut data = Vec::new();
        self.sequence(chromosome, start, stop, &mut data)?;
        Ok(data)
    }

    fn get_contig_list(&self) -> Vec<(String, u64)>;
    fn get_contig_length(&self, chromosome: &str) -> Option<u64>;
}

impl<R: Seek + Read + fmt::Debug> GenomeSequence for IndexedReader<R> {
    fn sequence(
        &mut self,
        chromosome: &str,
        start: u64,
        stop: u64,
        text: &mut Vec<u8>,
    ) -> Result<(), LiftOverError> {
        self.fetch(chromosome, start, stop)
            .map_err(|e| LiftOverError::UnknownSequenceError(chromosome.to_string(), e))?;
        //trace!("Sequence fetch: {}:{}-{}", chromosome, start, stop);
        self.read(text)?;
        //trace!("fetch OK");
        for x in text.iter_mut() {
            *x = match x {
                b'a' => b'A',
                b't' => b'T',
                b'c' => b'C',
                b'g' => b'G',
                _ => *x,
            };
        }
        Ok(())
    }

    fn get_contig_list(&self) -> Vec<(String, u64)> {
        self.index
            .sequences()
            .iter()
            .map(|x| (x.name.clone(), x.len))
            .collect()
    }

    fn get_contig_length(&self, chromosome: &str) -> Option<u64> {
        self.index
            .sequences()
            .iter()
            .find(|x| x.name == chromosome)
            .map(|x| x.len)
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct Variant {
    pub chromosome: String,
    /// zero based position
    pub position: u64,
    pub reference: Vec<u8>,
    pub alternative: Vec<Vec<u8>>,
}

impl Variant {
    pub fn new(
        chromosome: &str,
        position: u64,
        reference: &[u8],
        alternative: &[&[u8]],
    ) -> Variant {
        Variant {
            chromosome: chromosome.to_string(),
            position,
            reference: reference.to_vec(),
            alternative: alternative.iter().map(|x| x.to_vec()).collect(),
        }
    }
}

impl<R: VCFRecord> From<&R> for Variant {
    fn from(record: &R) -> Self {
        Variant {
            chromosome: str::from_utf8(record.contig()).unwrap().to_string(),
            position: record.position() - 1,
            reference: record.reference().to_vec(),
            alternative: if record.alternative() == [&b"."[..]] {
                vec![record.reference().to_vec()]
            } else {
                record.alternative().iter().map(|x| x.to_vec()).collect()
            },
        }
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut complement: Vec<_> = seq.iter().map(|x| reverse_acid(*x)).collect();

    complement.reverse();
    complement
}

pub fn reverse_acid(acid: u8) -> u8 {
    match acid {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'n' => b'n',
        b'N' => b'N',
        b'*' => b'*',
        _ => b'X',
    }
}

pub fn chromosome_priority(name: &str) -> u64 {
    let base_name = if name.starts_with("chr") {
        &name[3..]
    } else {
        name
    };
    match base_name {
        "X" => 23,
        "Y" => 24,
        "M" => 25,
        _ => {
            if let Ok(chrom_num) = base_name.parse::<u64>() {
                chrom_num
            } else {
                26
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_chromosome_priority() {
        assert_eq!(chromosome_priority("chr1"), 1);
        assert_eq!(chromosome_priority("1"), 1);
        assert_eq!(chromosome_priority("chr10"), 10);
        assert_eq!(chromosome_priority("10"), 10);
        assert_eq!(chromosome_priority("chrX"), 23);
        assert_eq!(chromosome_priority("X"), 23);
        assert_eq!(chromosome_priority("chrY"), 24);
        assert_eq!(chromosome_priority("Y"), 24);
        assert_eq!(chromosome_priority("chrM"), 25);
        assert_eq!(chromosome_priority("M"), 25);
        assert_eq!(chromosome_priority("hoge"), 26);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b""), b"");
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"T"), b"A");
        assert_eq!(reverse_complement(b"C"), b"G");
        assert_eq!(reverse_complement(b"G"), b"C");
        assert_eq!(reverse_complement(b"AT"), b"AT");
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
        assert_eq!(reverse_complement(b"ABCGC"), b"GCGXT");
    }
}
