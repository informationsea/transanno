use super::{Feature, FeatureType, Gene, GeneParseError, GeneStrand, GroupedReader, Transcript};
use indexmap::IndexMap;
use log::{error, trace};
use nom::branch::alt;
use nom::bytes::complete::{is_not, tag};
use nom::character::complete::line_ending;
use nom::multi::separated_list;
use nom::sequence::{delimited, separated_pair};
use once_cell::sync::Lazy;
use regex::Regex;
use std::collections::HashMap;
use std::fmt;
use std::io;
use std::iter::Iterator;
use std::str::FromStr;
use std::string::ToString;

static DIGITS: Lazy<Regex> = Lazy::new(|| Regex::new("^\\d+$").unwrap());

#[derive(Debug)]
pub struct GtfReader<R: io::BufRead> {
    reader: R,
    line_number: u64,
}

impl<R: io::BufRead> GtfReader<R> {
    pub fn new(reader: R) -> GtfReader<R> {
        GtfReader {
            reader,
            line_number: 0,
        }
    }
}

impl<R: io::BufRead> Iterator for GtfReader<R> {
    type Item = Result<GtfRecord, GeneParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let mut line = String::new();
            self.line_number += 1;

            return match self.reader.read_line(&mut line) {
                Ok(0) => None,
                Ok(_) => {
                    if line.starts_with('#') {
                        continue;
                    } else {
                        match line.parse::<GtfRecord>() {
                            Ok(feature) => Some(Ok(feature)),
                            Err(e) => Some(Err(GeneParseError::ParseErrorAtLine(
                                self.line_number,
                                Box::new(e),
                            ))),
                        }
                    }
                }
                Err(e) => Some(Err(e.into())),
            };
        }
    }
}

pub struct GtfGroupedReader<R: io::BufRead> {
    reader: GtfReader<R>,
    last_item: Option<GtfRecord>,
}

impl<R: io::BufRead> GtfGroupedReader<R> {
    pub fn new(reader: GtfReader<R>) -> GtfGroupedReader<R> {
        GtfGroupedReader {
            reader,
            last_item: None,
        }
    }
}

impl<R: io::BufRead> Iterator for GtfGroupedReader<R> {
    type Item = Result<Gene<GtfRecord, GtfRecord, GtfRecord>, GeneParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut records = Vec::new();
        if let Some(v) = self.last_item.take() {
            records.push(v);
        }

        while let Some(v) = self.reader.next() {
            match v {
                Ok(v) => {
                    if v.record_type == "gene" && !records.is_empty() {
                        // must be gene
                        self.last_item = Some(v);
                        break;
                    } else {
                        records.push(v);
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        if records.is_empty() {
            return None;
        }

        Some(group_gtf(records))
    }
}

impl<R: io::BufRead> GroupedReader<GtfRecord, GtfRecord, GtfRecord> for GtfGroupedReader<R> {}

fn group_gtf(
    mut records: Vec<GtfRecord>,
) -> Result<Gene<GtfRecord, GtfRecord, GtfRecord>, GeneParseError> {
    let gene_candidates: Vec<_> = records
        .iter()
        .enumerate()
        .filter(|x| x.1.record_type == "gene")
        .collect();

    if gene_candidates.len() != 1 {
        error!(
            "Cannot group by gene: multiple genes or no gene in a group: {:?}",
            gene_candidates
        );
        return Err(GeneParseError::ParseError.into());
    }

    let gene_index = gene_candidates[0].0;

    let mut gene = Gene {
        original_record: records.remove(gene_index),
        transcripts: Vec::new(),
    };
    trace!("grouping gene: {:?}", gene);

    if !gene.original_record.attributes.contains_key("gene_id") {
        println!(
            "No ID was found for gene candidate record: {:?}",
            gene.original_record
        );
        return Err(GeneParseError::ParseError.into());
    }

    let transcript_indexes: Vec<usize> = records
        .iter()
        .enumerate()
        .filter(|x| x.1.attributes.contains_key("transcript_id"))
        .filter(|x| x.1.record_type == "transcript")
        .filter(|x| x.1.attributes.get("gene_id") == gene.original_record.attributes.get("gene_id"))
        .map(|x| x.0)
        .rev()
        .collect();
    let mut transcripts: HashMap<_, _> = transcript_indexes
        .iter()
        .map(|x| records.remove(*x))
        .map(|x| {
            (
                x.attributes.get("transcript_id").unwrap().clone(),
                Transcript {
                    original_record: x,
                    children: Vec::new(),
                },
            )
        })
        .collect();
    //println!("transcripts: {:?}", transcripts.keys());

    for one_record in records {
        let parent_id = one_record
            .attributes
            .get("transcript_id")
            .expect("No transcript ID");
        if let Some(parent) = transcripts.get_mut(parent_id) {
            parent.children.push(one_record);
        } else {
            error!("no parent found: {:?}", parent_id);
            return Err(GeneParseError::ParseError.into());
        }
    }

    for (_, v) in transcripts {
        gene.transcripts.push(v);
    }

    Ok(gene)
}

#[derive(Debug, Clone, PartialEq)]
pub struct GtfRecord {
    pub seq_id: String,
    pub source: String,
    pub feature_type: FeatureType,
    pub record_type: String,
    pub start: u64,
    pub end: u64,
    pub score: String,
    pub strand: GeneStrand,
    pub phase: Option<u8>,
    pub attributes: IndexMap<String, Vec<String>>,
}

impl Feature for GtfRecord {
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

    fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .get(key)
            .and_then(|x| x.get(0).map(|y| y.as_str()))
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

    fn set_attribute(&mut self, key: &str, value: &str) {
        self.attributes
            .insert(key.to_string(), vec![value.to_string()]);
    }
}

impl fmt::Display for GtfRecord {
    fn fmt(&self, writer: &mut fmt::Formatter) -> fmt::Result {
        write!(writer, "{}\t", self.seq_id)?;
        write!(writer, "{}\t", self.source)?;
        write!(writer, "{}\t", self.record_type)?;
        write!(writer, "{}\t", self.start)?;
        write!(writer, "{}\t", self.end)?;
        write!(writer, "{}\t", self.score)?;
        write!(
            writer,
            "{}\t",
            match self.strand {
                GeneStrand::Forward => "+",
                GeneStrand::Reverse => "-",
                _ => ".",
            }
        )?;
        if let Some(phase) = self.phase {
            write!(writer, "{}\t", phase)?;
        } else {
            write!(writer, ".\t")?;
        }

        if self.attributes.is_empty() {
            write!(writer, ".")?;
        } else {
            for (i, (k, v)) in self.attributes.iter().enumerate() {
                for (j, item) in v.iter().enumerate() {
                    if i != 0 || j != 0 {
                        write!(writer, " ")?;
                    }
                    if DIGITS.is_match(item) {
                        write!(writer, "{} {};", k, item)?;
                    } else {
                        write!(writer, "{} \"{}\";", k, item)?;
                    }
                }
            }
        }
        writeln!(writer)?;

        Ok(())
    }
}

impl FromStr for GtfRecord {
    type Err = GeneParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        let (input, seq_id) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, source) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, record_type) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, start) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, end) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, score) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, strand) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, phase) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, attributes) = separated_list(
            tag("; "),
            separated_pair(
                is_not(" \t\r\n\";"),
                tag(" "),
                alt((
                    delimited(tag("\""), is_not("\""), tag("\"")),
                    is_not(" \t\r\n\";"),
                )),
            ),
        )(input)?;
        let (input, _) = tag(";")(input)?;
        line_ending(input)?;

        Ok(GtfRecord {
            seq_id: seq_id.to_string(),
            source: source.to_string(),
            feature_type: record_type.parse()?,
            record_type: record_type.to_string(),
            start: start.parse()?,
            end: end.parse()?,
            score: score.to_string(),
            strand: match strand {
                "+" => GeneStrand::Forward,
                "-" => GeneStrand::Reverse,
                _ => GeneStrand::Unknown,
            },
            phase: match phase {
                "0" => Some(0),
                "1" => Some(1),
                "2" => Some(2),
                _ => None,
            },
            attributes: attributes.into_iter().fold(
                IndexMap::new(),
                |mut acc: IndexMap<_, _>, (k, v)| {
                    if let Some(l) = acc.get_mut(k) {
                        l.push(v.to_string());
                    } else {
                        acc.insert(k.to_string(), vec![v.to_string()]);
                    }
                    acc
                },
            ),
        })
    }
}

#[cfg(test)]
#[allow(clippy::unreadable_literal)]
mod test {
    use super::*;
    use indexmap::indexmap;
    use lazy_static::lazy_static;

    static CHEK2P4_GENE_LINE: &str = r#"chr22	HAVANA	gene	16503304	16509475	.	+	.	gene_id "ENSG00000235969.1"; gene_type "unprocessed_pseudogene"; gene_name "CHEK2P4"; level 2; hgnc_id "HGNC:43580"; havana_gene "OTTHUMG00000140400.1";
"#;

    static CHEK2_CDS4_LINE: &str = r#"chr22	HAVANA	CDS	28699838	28699937	.	-	1	gene_id "ENSG00000183765.22"; transcript_id "ENST00000404276.6"; gene_type "protein_coding"; gene_name "CHEK2"; transcript_type "protein_coding"; transcript_name "CHEK2-207"; exon_number 9; exon_id "ENSE00003694540.1"; level 2; protein_id "ENSP00000385747.1"; transcript_support_level "1"; hgnc_id "HGNC:16627"; tag "CAGE_supported_TSS"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS13843.1"; havana_gene "OTTHUMG00000151023.24"; havana_transcript "OTTHUMT00000500899.1";
"#;

    lazy_static! {
        static ref CHEK2P4_GENE: GtfRecord = GtfRecord {
            seq_id: "chr22".to_string(),
            source: "HAVANA".to_string(),
            feature_type: FeatureType::Gene,
            record_type: "gene".to_string(),
            start: 16503304,
            end: 16509475,
            score: ".".to_string(),
            strand: GeneStrand::Forward,
            phase: None,
            attributes: indexmap! {
                "gene_id".to_string() => vec!["ENSG00000235969.1".to_string()],
                "gene_type".to_string() => vec!["unprocessed_pseudogene".to_string()],
                "gene_name".to_string() => vec!["CHEK2P4".to_string()],
                "level".to_string() => vec!["2".to_string()],
                "hgnc_id".to_string() => vec!["HGNC:43580".to_string()],
                "havana_gene".to_string() => vec!["OTTHUMG00000140400.1".to_string()],
            },
        };
        static ref CHEK2_CDS4: GtfRecord = GtfRecord {
            seq_id: "chr22".to_string(),
            source: "HAVANA".to_string(),
            feature_type: FeatureType::CDS,
            record_type: "CDS".to_string(),
            start: 28699838,
            end: 28699937,
            score: ".".to_string(),
            strand: GeneStrand::Reverse,
            phase: Some(1),
            attributes: indexmap! {
                "gene_id".to_string() => vec!["ENSG00000183765.22".to_string()],
                "transcript_id".to_string() => vec!["ENST00000404276.6".to_string()],
                "gene_type".to_string() => vec!["protein_coding".to_string()],
                "gene_name".to_string() => vec!["CHEK2".to_string()],
                "transcript_type".to_string() => vec!["protein_coding".to_string()],
                "transcript_name".to_string() => vec!["CHEK2-207".to_string()],
                "exon_number".to_string() => vec!["9".to_string()],
                "exon_id".to_string() => vec!["ENSE00003694540.1".to_string()],
                "level".to_string() => vec!["2".to_string()],
                "protein_id".to_string() => vec!["ENSP00000385747.1".to_string()],
                "transcript_support_level".to_string() => vec!["1".to_string()],
                "hgnc_id".to_string() => vec!["HGNC:16627".to_string()],
                "tag".to_string() => vec!["CAGE_supported_TSS".to_string(), "basic".to_string(), "MANE_Select".to_string(), "appris_principal_3".to_string(), "CCDS".to_string()],
                "ccdsid".to_string() => vec!["CCDS13843.1".to_string()],
                "havana_gene".to_string() => vec!["OTTHUMG00000151023.24".to_string()],
                "havana_transcript".to_string() => vec!["OTTHUMT00000500899.1".to_string()],
            },
        };
    }

    #[allow(clippy::unreadable_literal)]
    #[test]
    fn test_gtf_parse() {
        assert_eq!(
            CHEK2P4_GENE_LINE.parse::<GtfRecord>().unwrap(),
            *CHEK2P4_GENE
        );

        assert_eq!(CHEK2_CDS4_LINE.parse::<GtfRecord>().unwrap(), *CHEK2_CDS4);
    }

    #[allow(clippy::unreadable_literal)]
    #[test]
    fn test_gtf_to_string() {
        assert_eq!(CHEK2P4_GENE_LINE, CHEK2P4_GENE.to_string());
        //assert_eq!(CHEK2_CDS4_LINE, CHEK2_CDS4.to_string());
    }

    #[test]
    fn test_gtf_reader() {
        let test_data =
            autocompress::Decoder::suggest(
                &include_bytes!(
                    "../../testfiles/GENCODE/gencode.v33.basic.annotation.chr22.gtf.xz"
                )[..],
            )
            .unwrap();
        let test_reader = io::BufReader::new(test_data);
        let gtf_reader = GtfReader::new(test_reader);
        let values: Vec<_> = gtf_reader.map(|x| x.unwrap()).collect();
        assert_eq!(values[352], *CHEK2P4_GENE);
        assert_eq!(values[13003], *CHEK2_CDS4);
        assert_eq!(values.len(), 37907);
    }

    #[test]
    fn test_gtf_grouped_reader() {
        let test_data =
            autocompress::Decoder::suggest(
                &include_bytes!(
                    "../../testfiles/GENCODE/gencode.v33.basic.annotation.chr22.gtf.xz"
                )[..],
            )
            .unwrap();
        let test_reader = io::BufReader::new(test_data);
        let gtf_reader = GtfReader::new(test_reader);
        let grouped_gtf_reader = GtfGroupedReader::new(gtf_reader);
        let values: Vec<_> = grouped_gtf_reader.collect();

        assert_eq!(values[48].as_ref().unwrap().original_record, *CHEK2P4_GENE);
        assert_eq!(values.len(), 1388);
        // TODO: Check more details
    }
}
