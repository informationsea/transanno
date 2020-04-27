mod error;
pub mod gff3;
pub mod gtf;
pub mod refgene;

pub use error::{GeneParseError, GeneParseErrorKind};

use std::fmt::{self, Display};
use std::iter::Iterator;
use std::str::FromStr;

pub trait GroupedReader<G: Feature, T: Feature, F: Feature>:
    Iterator<Item = Result<Gene<G, T, F>, GeneParseError>>
{
}

#[derive(Clone, Copy, Hash, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum FeatureType {
    Gene,
    Transcript,
    Exon,
    CDS,
    ThreePrimeUTR,
    FivePrimeUTR,
    StartCodon,
    StopCodon,
    UTR,
    Other,
}

impl FromStr for FeatureType {
    type Err = GeneParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "transcript" | "pseudogenic_transcript" => FeatureType::Transcript,
            "gene" | "pseudogene" => FeatureType::Gene,
            "exon" => FeatureType::Exon,
            "CDS" => FeatureType::CDS,
            "three_prime_UTR" => FeatureType::ThreePrimeUTR,
            "five_prime_UTR" => FeatureType::FivePrimeUTR,
            "start_codon" => FeatureType::StartCodon,
            "stop_codon" => FeatureType::StopCodon,
            "UTR" => FeatureType::UTR,
            _ => FeatureType::Other,
        })
    }
}

#[derive(Clone, Copy, Hash, Debug, PartialEq, PartialOrd, Eq, Ord)]
pub enum GeneStrand {
    Forward,
    Reverse,
    Unknown,
}

impl std::ops::Mul for GeneStrand {
    type Output = GeneStrand;
    fn mul(self, rhs: Self) -> Self::Output {
        match self {
            GeneStrand::Forward => match rhs {
                GeneStrand::Forward => GeneStrand::Forward,
                GeneStrand::Reverse => GeneStrand::Reverse,
                GeneStrand::Unknown => GeneStrand::Unknown,
            },
            GeneStrand::Reverse => match rhs {
                GeneStrand::Forward => GeneStrand::Reverse,
                GeneStrand::Reverse => GeneStrand::Forward,
                GeneStrand::Unknown => GeneStrand::Unknown,
            },
            GeneStrand::Unknown => GeneStrand::Unknown,
        }
    }
}

pub trait Feature: fmt::Debug + Clone + Sized + PartialEq {
    /// Chromosome name    
    fn seq_id(&self) -> &str;
    /// Feature type
    fn feature_type(&self) -> FeatureType;
    /// 1-based start position (included)
    fn start(&self) -> u64;
    /// 1-based end position (included)
    fn end(&self) -> u64;
    /// strand
    fn strand(&self) -> GeneStrand;
    /// attributes
    fn attribute(&self, key: &str) -> Option<&str>;

    /// 0-based range
    fn range(&self) -> std::ops::Range<u64> {
        (self.start() - 1)..self.end()
    }

    fn len(&self) -> u64 {
        self.end() + 1 - self.start()
    }

    fn is_empty(&self) -> bool {
        self.start() == (self.end() + 1)
    }

    fn seq_id_mut(&mut self) -> &mut String;
    fn start_mut(&mut self) -> &mut u64;
    fn end_mut(&mut self) -> &mut u64;
    fn strand_mut(&mut self) -> &mut GeneStrand;
    fn set_attribute(&mut self, key: &str, value: &str);

    fn set_range(&mut self, range: std::ops::Range<u64>) {
        *self.start_mut() = range.start + 1;
        *self.end_mut() = range.end;
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Gene<G: Feature, T: Feature, F: Feature> {
    pub original_record: G,
    pub transcripts: Vec<Transcript<T, F>>,
}

impl<G: Feature + Display, T: Feature + Display, F: Feature + Display> Display for Gene<G, T, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.original_record)?;
        for one in self.transcripts.iter() {
            write!(f, "{}", one)?;
        }

        Ok(())
    }
}

impl<G: Feature, T: Feature, F: Feature> Feature for Gene<G, T, F> {
    fn seq_id(&self) -> &str {
        self.original_record.seq_id()
    }

    fn feature_type(&self) -> FeatureType {
        self.original_record.feature_type()
    }
    fn start(&self) -> u64 {
        self.original_record.start()
    }

    fn end(&self) -> u64 {
        self.original_record.end()
    }

    fn strand(&self) -> GeneStrand {
        self.original_record.strand()
    }

    fn attribute(&self, key: &str) -> Option<&str> {
        self.original_record.attribute(key)
    }

    fn seq_id_mut(&mut self) -> &mut String {
        self.original_record.seq_id_mut()
    }

    fn start_mut(&mut self) -> &mut u64 {
        self.original_record.start_mut()
    }

    fn end_mut(&mut self) -> &mut u64 {
        self.original_record.end_mut()
    }

    fn strand_mut(&mut self) -> &mut GeneStrand {
        self.original_record.strand_mut()
    }

    fn set_attribute(&mut self, key: &str, value: &str) {
        self.original_record.set_attribute(key, value)
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Transcript<T: Feature, F: Feature> {
    pub original_record: T,
    pub children: Vec<F>,
}

impl<T: Feature + Display, F: Feature + Display> Display for Transcript<T, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.original_record)?;
        for one in self.children.iter() {
            write!(f, "{}", one)?;
        }
        Ok(())
    }
}

impl<T: Feature, F: Feature> Feature for Transcript<T, F> {
    fn seq_id(&self) -> &str {
        self.original_record.seq_id()
    }

    fn feature_type(&self) -> FeatureType {
        self.original_record.feature_type()
    }

    fn start(&self) -> u64 {
        self.original_record.start()
    }

    fn end(&self) -> u64 {
        self.original_record.end()
    }

    fn strand(&self) -> GeneStrand {
        self.original_record.strand()
    }

    fn attribute(&self, key: &str) -> Option<&str> {
        self.original_record.attribute(key)
    }

    fn seq_id_mut(&mut self) -> &mut String {
        self.original_record.seq_id_mut()
    }

    fn start_mut(&mut self) -> &mut u64 {
        self.original_record.start_mut()
    }

    fn end_mut(&mut self) -> &mut u64 {
        self.original_record.end_mut()
    }

    fn strand_mut(&mut self) -> &mut GeneStrand {
        self.original_record.strand_mut()
    }

    fn set_attribute(&mut self, key: &str, value: &str) {
        self.original_record.set_attribute(key, value)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[derive(Debug, Clone, PartialEq, Eq, PartialOrd)]
    struct TestFeature {
        seq_id: String,
        feature_type: FeatureType,
        start: u64,
        end: u64,
        strand: GeneStrand,
    }

    impl Feature for TestFeature {
        fn seq_id(&self) -> &str {
            &self.seq_id
        }

        fn feature_type(&self) -> FeatureType {
            self.feature_type
        }

        fn start(&self) -> u64 {
            self.start
        }

        fn end(&self) -> u64 {
            self.end
        }

        fn strand(&self) -> GeneStrand {
            self.strand
        }
        fn attribute(&self, _key: &str) -> Option<&str> {
            unimplemented!()
        }

        fn seq_id_mut(&mut self) -> &mut String {
            &mut self.seq_id
        }
        fn start_mut(&mut self) -> &mut u64 {
            &mut self.start
        }
        fn end_mut(&mut self) -> &mut u64 {
            &mut self.end
        }
        fn strand_mut(&mut self) -> &mut GeneStrand {
            &mut self.strand
        }
        fn set_attribute(&mut self, _key: &str, _value: &str) {
            unimplemented!()
        }
    }

    #[test]
    fn test_feature() {
        let mut x = TestFeature {
            seq_id: "chr22".to_string(),
            start: 100,
            end: 120,
            strand: GeneStrand::Forward,
            feature_type: FeatureType::Exon,
        };

        assert_eq!(x.len(), 21);
        assert_eq!(x.range(), 99..120);
        assert_eq!(x.is_empty(), false);

        x.set_range(10..10);
        assert_eq!(
            x,
            TestFeature {
                seq_id: "chr22".to_string(),
                start: 11,
                end: 10,
                strand: GeneStrand::Forward,
                feature_type: FeatureType::Exon,
            }
        );
        assert_eq!(x.is_empty(), true);
        assert_eq!(x.len(), 0);
    }
}
