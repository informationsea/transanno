use super::{GeneParseError, GeneStrand};
use std::fmt::{self, Display};
use std::io::BufRead;
use std::str::FromStr;

#[derive(Debug)]
pub struct RefGeneReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> RefGeneReader<R> {
    pub fn new(reader: R) -> Self {
        RefGeneReader { reader }
    }
}

impl<R: BufRead> Iterator for RefGeneReader<R> {
    type Item = Result<RefGeneEntry, GeneParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => return None,
            Err(e) => return Some(Err(e.into())),
            _ => (),
        }

        Some(line.parse())
    }
}

#[derive(Debug, Clone, PartialEq, PartialOrd, Eq, Ord)]
pub enum CDSStat {
    None,
    Unk,
    Incmpl,
    Cmpl,
}

impl Display for CDSStat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CDSStat::None => write!(f, "none"),
            CDSStat::Unk => write!(f, "unk"),
            CDSStat::Incmpl => write!(f, "incmpl"),
            CDSStat::Cmpl => write!(f, "cmpl"),
        }
    }
}

impl FromStr for CDSStat {
    type Err = GeneParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "none" => Ok(CDSStat::None),
            "unk" => Ok(CDSStat::Unk),
            "incmpl" => Ok(CDSStat::Incmpl),
            "cmpl" => Ok(CDSStat::Cmpl),
            _ => Err(GeneParseError::CDSStatParseError),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RefGeneEntry {
    pub bin: u32,
    pub name: String,
    pub chrom: String,
    pub strand: GeneStrand,
    pub tx_start: u64,
    pub tx_end: u64,
    pub cds_start: u64,
    pub cds_end: u64,
    pub exon_count: u32,
    pub exon_starts: Vec<u64>,
    pub exon_ends: Vec<u64>,
    pub score: u64,
    pub name2: String,
    pub cds_start_stat: CDSStat,
    pub cds_end_stat: CDSStat,
    pub exon_frames: Vec<i32>,
}

impl FromStr for RefGeneEntry {
    type Err = GeneParseError;
    fn from_str(line: &str) -> Result<Self, Self::Err> {
        let elements: Vec<_> = line.trim_end().split('\t').collect();
        if elements.len() != 16 {
            return Err(GeneParseError::RefGeneInvalidColumnNumber);
        }

        Ok(RefGeneEntry {
            bin: elements[0].parse()?,
            name: elements[1].to_string(),
            chrom: elements[2].to_string(),
            strand: match elements[3] {
                "+" => GeneStrand::Forward,
                "-" => GeneStrand::Reverse,
                _ => GeneStrand::Unknown,
            },
            tx_start: elements[4].parse()?,
            tx_end: elements[5].parse()?,
            cds_start: elements[6].parse()?,
            cds_end: elements[7].parse()?,
            exon_count: elements[8].parse()?,
            exon_starts: elements[9].split(',').filter(|x| !x.is_empty()).try_fold(
                Vec::new(),
                |mut acc, x| -> Result<Vec<u64>, GeneParseError> {
                    acc.push(x.parse()?);
                    Ok(acc)
                },
            )?,
            exon_ends: elements[10].split(',').filter(|x| !x.is_empty()).try_fold(
                Vec::new(),
                |mut acc, x| -> Result<Vec<u64>, GeneParseError> {
                    acc.push(x.parse()?);
                    Ok(acc)
                },
            )?,
            score: elements[11].parse()?,
            name2: elements[12].to_string(),
            cds_start_stat: elements[13].parse()?,
            cds_end_stat: elements[14].parse()?,
            exon_frames: elements[15].split(',').filter(|x| !x.is_empty()).try_fold(
                Vec::new(),
                |mut acc, x| -> Result<Vec<i32>, GeneParseError> {
                    acc.push(x.parse()?);
                    Ok(acc)
                },
            )?,
        })
    }
}

// TODO: Implement Grouped RefGene Reader

// #[derive(Debug, Clone)]
// pub struct RefGeneTranscript {
//     pub seq_id: String,
//     pub source: String,
//     pub start: u64,
//     pub end: u64,
//     pub strand: Option<GeneStrand>,
// }

// impl Feature for RefGeneTranscript {
//     fn seq_id(&self) -> &str {
//         &self.seq_id
//     }

//     fn record_type(&self) -> &str {
//         "transcript"
//     }

//     fn start(&self) -> u64 {
//         self.start
//     }

//     fn end(&self) -> u64 {
//         self.end
//     }

//     fn strand(&self) -> Option<GeneStrand> {
//         self.strand
//     }

//     fn seq_id_mut(&mut self) -> &mut String {
//         &mut self.seq_id
//     }

//     fn start_mut(&mut self) -> &mut u64 {
//         &mut self.start
//     }

//     fn end_mut(&mut self) -> &mut u64 {
//         &mut self.end
//     }

//     fn strand_mut(&mut self) -> &mut Option<GeneStrand> {
//         &mut self.strand
//     }
// }

// #[derive(Debug, Clone)]
// pub struct RefGeneFeature {
//     pub seq_id: String,
//     pub source: String,
//     pub record_type: String,
//     pub start: u64,
//     pub end: u64,
//     pub score: String,
//     pub strand: Option<GeneStrand>,
//     pub phase: Option<u8>,
// }

// impl Feature for RefGeneFeature {
//     fn seq_id(&self) -> &str {
//         &self.seq_id
//     }

//     fn record_type(&self) -> &str {
//         &self.record_type
//     }

//     fn start(&self) -> u64 {
//         self.start
//     }

//     fn end(&self) -> u64 {
//         self.end
//     }

//     fn strand(&self) -> Option<GeneStrand> {
//         self.strand
//     }

//     fn seq_id_mut(&mut self) -> &mut String {
//         &mut self.seq_id
//     }

//     fn start_mut(&mut self) -> &mut u64 {
//         &mut self.start
//     }

//     fn end_mut(&mut self) -> &mut u64 {
//         &mut self.end
//     }

//     fn strand_mut(&mut self) -> &mut Option<GeneStrand> {
//         &mut self.strand
//     }
// }

#[allow(clippy::unreadable_literal)]
#[cfg(test)]
mod test {
    use once_cell::sync::Lazy;

    use super::*;

    static CHEK2_LINE: &str = r#"100	NM_145862	chr22	-	28687742	28741834	28687896	28734721	14	28687742,28689134,28694031,28695126,28695709,28699837,28703504,28710005,28711908,28719394,28724976,28725242,28734402,28741768,	28687986,28689215,28694117,28695242,28695873,28699937,28703566,28710059,28712017,28719485,28725124,28725367,28734727,28741834,	0	CHEK2	cmpl	cmpl	0,0,1,2,0,2,0,0,2,1,0,1,0,-1,
"#;

    static CHEK2: Lazy<RefGeneEntry> = Lazy::new(|| RefGeneEntry {
        bin: 100,
        name: "NM_145862".to_string(),
        chrom: "chr22".to_string(),
        strand: GeneStrand::Reverse,
        tx_start: 28687742,
        tx_end: 28741834,
        cds_start: 28687896,
        cds_end: 28734721,
        exon_count: 14,
        exon_starts: vec![
            28687742, 28689134, 28694031, 28695126, 28695709, 28699837, 28703504, 28710005,
            28711908, 28719394, 28724976, 28725242, 28734402, 28741768,
        ],
        exon_ends: vec![
            28687986, 28689215, 28694117, 28695242, 28695873, 28699937, 28703566, 28710059,
            28712017, 28719485, 28725124, 28725367, 28734727, 28741834,
        ],
        score: 0,
        name2: "CHEK2".to_string(),
        cds_start_stat: CDSStat::Cmpl,
        cds_end_stat: CDSStat::Cmpl,
        exon_frames: vec![0, 0, 1, 2, 0, 2, 0, 0, 2, 1, 0, 1, 0, -1],
    });

    #[test]
    fn test_refgene_reader() -> Result<(), GeneParseError> {
        let parsed: RefGeneEntry = CHEK2_LINE.parse()?;
        assert_eq!(parsed, *CHEK2);
        Ok(())
    }
}
