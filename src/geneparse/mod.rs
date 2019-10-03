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

#[derive(Clone, Copy, Hash, Debug, PartialEq, PartialOrd)]
pub enum FeatureType {
    Transcript,
    Gene,
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

pub trait Feature: fmt::Debug + Clone + Sized + PartialEq + Display {
    fn seq_id(&self) -> &str;
    fn feature_type(&self) -> FeatureType;
    fn start(&self) -> u64;
    fn end(&self) -> u64;
    fn strand(&self) -> GeneStrand;
    fn attribute(&self, key: &str) -> Option<&str>;

    fn seq_id_mut(&mut self) -> &mut String;
    fn start_mut(&mut self) -> &mut u64;
    fn end_mut(&mut self) -> &mut u64;
    fn strand_mut(&mut self) -> &mut GeneStrand;
    fn set_attribute(&mut self, key: &str, value: &str);
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Gene<G: Feature, T: Feature, F: Feature> {
    pub original_record: G,
    pub transcripts: Vec<Transcript<T, F>>,
}

impl<G: Feature, T: Feature, F: Feature> Display for Gene<G, T, F> {
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

impl<T: Feature, F: Feature> Display for Transcript<T, F> {
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
