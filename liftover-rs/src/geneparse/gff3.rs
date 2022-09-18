use super::{Feature, FeatureType, Gene, GeneParseError, GeneStrand, GroupedReader, Transcript};
use indexmap::IndexMap;
use log::{error, trace};
use nom::bytes::complete::{is_not, tag};
use nom::character::complete::line_ending;
use nom::multi::separated_list;
use nom::sequence::separated_pair;
use std::collections::HashMap;
use std::fmt;
use std::io;
use std::iter::Iterator;
use std::str::FromStr;

#[derive(Debug)]
pub struct Gff3Reader<R: io::BufRead> {
    reader: R,
    line_number: u64,
}

impl<R: io::BufRead> Gff3Reader<R> {
    pub fn new(reader: R) -> Gff3Reader<R> {
        Gff3Reader {
            reader,
            line_number: 0,
        }
    }
}

impl<R: io::BufRead> Iterator for Gff3Reader<R> {
    type Item = Result<Gff3Record, GeneParseError>;

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
                        match line.parse::<Gff3Record>() {
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

pub struct Gff3GroupedReader<R: io::BufRead> {
    reader: Gff3Reader<R>,
    last_item: Option<Gff3Record>,
}

impl<R: io::BufRead> Gff3GroupedReader<R> {
    pub fn new(reader: Gff3Reader<R>) -> Gff3GroupedReader<R> {
        Gff3GroupedReader {
            reader,
            last_item: None,
        }
    }
}

impl<R: io::BufRead> GroupedReader<Gff3Record, Gff3Record, Gff3Record> for Gff3GroupedReader<R> {}

impl<R: io::BufRead> Iterator for Gff3GroupedReader<R> {
    type Item = Result<Gene<Gff3Record, Gff3Record, Gff3Record>, GeneParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut records = Vec::new();
        if let Some(v) = self.last_item.take() {
            records.push(v);
        }

        while let Some(v) = self.reader.next() {
            match v {
                Ok(v) => {
                    let parent = v.attributes.get("Parent");
                    let id = v.attributes.get("ID");
                    if parent.is_none() && id.is_none() {
                        // skip "biological_region"
                        continue;
                    } else if parent.is_none() && !records.is_empty() {
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

        Some(group_gff3(records))
    }
}

fn group_gff3(
    mut records: Vec<Gff3Record>,
) -> Result<Gene<Gff3Record, Gff3Record, Gff3Record>, GeneParseError> {
    let gene_candidates: Vec<_> = records
        .iter()
        .enumerate()
        .filter(|x| !x.1.attributes.contains_key("Parent"))
        .collect();

    if gene_candidates.len() != 1 {
        error!(
            "Cannot group by gene: multiple genes or no gene in a group: {:?}",
            gene_candidates
        );
        return Err(GeneParseError::ParseError);
    }

    let gene_index = gene_candidates[0].0;

    let mut gene = Gene {
        original_record: records.remove(gene_index),
        transcripts: Vec::new(),
    };
    trace!("grouping gene: {:?}", gene);

    if !gene.original_record.attributes.contains_key("ID") {
        error!(
            "No ID was found for gene candidate record: {:?}",
            gene.original_record
        );
        return Err(GeneParseError::ParseError);
    }

    let transcript_indexes: Vec<usize> = records
        .iter()
        .enumerate()
        .filter(|x| x.1.attributes.contains_key("ID"))
        .filter(|x| x.1.attributes.get("Parent") == gene.original_record.attributes.get("ID"))
        .map(|x| x.0)
        .rev()
        .collect();
    let mut transcripts: HashMap<_, _> = transcript_indexes
        .iter()
        .map(|x| records.remove(*x))
        .map(|x| {
            (
                x.attributes.get("ID").unwrap().clone(),
                Transcript {
                    original_record: x,
                    children: Vec::new(),
                },
            )
        })
        .collect();
    trace!("transcripts: {:?}", transcripts.keys());

    let mut parent_remap = HashMap::new();

    for one_record in records {
        let parent_id = one_record.attributes.get("Parent").expect("No parent ID");
        if let Some(parent) = transcripts.get_mut(parent_id) {
            if let Some(id) = one_record.attributes.get("ID") {
                parent_remap.insert(id.to_string(), parent_id.to_string());
            }

            parent.children.push(one_record);
        } else if let Some(grand_parent_id) = parent_remap.get(parent_id) {
            if let Some(parent) = transcripts.get_mut(grand_parent_id) {
                parent.children.push(one_record);
            } else {
                println!("no parent found: {:?} {:?}", parent_id, one_record);
                // Skip no parent record (e.g. stop_codon_redefined_as_selenocysteine)
                //return Err(GeneParseErrorKind::ParseError.into());
            }
        } else {
            println!("no parent found: {:?} {:?}", parent_id, one_record);
            // Skip no parent record (e.g. stop_codon_redefined_as_selenocysteine)
            //return Err(GeneParseErrorKind::ParseError.into());
        }
    }

    for (_, v) in transcripts {
        gene.transcripts.push(v);
    }

    Ok(gene)
}

#[derive(Debug, Clone, PartialEq)]
pub struct Gff3Record {
    pub seq_id: String,
    pub source: String,
    pub feature_type: FeatureType,
    pub record_type: String,
    pub start: u64,
    pub end: u64,
    pub score: String,
    pub strand: GeneStrand,
    pub phase: Option<u8>,
    pub attributes: IndexMap<String, String>,
}

impl Feature for Gff3Record {
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
        self.attributes.get(key).map(|x| x.as_str())
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
        self.attributes.insert(key.to_string(), value.to_string());
    }
}

impl fmt::Display for Gff3Record {
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
                if i != 0 {
                    write!(writer, ";")?;
                }
                write!(writer, "{}={}", k, v)?;
            }
        }
        writeln!(writer)?;

        Ok(())
    }
}

impl FromStr for Gff3Record {
    type Err = GeneParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        let (input, seq_id) = is_not(" \t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, source) = is_not("\t\r\n")(input)?;
        let (input, _) = tag("\t")(input)?;
        let (input, record_type) = is_not("\t\r\n")(input)?;
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
            tag(";"),
            separated_pair(is_not("\r\n=;"), tag("="), is_not("\r\n=;")),
        )(input)?;
        line_ending(input)?;

        Ok(Gff3Record {
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
            attributes: attributes
                .iter()
                .map(|(k, v)| ((*k).to_string(), (*v).to_string()))
                .collect(),
        })
    }
}

#[allow(clippy::unreadable_literal)]
#[cfg(test)]
mod test {
    use super::*;
    use indexmap::indexmap;
    use lazy_static::lazy_static;

    lazy_static! {
        static ref CHEK2P4_GENE: Gff3Record = Gff3Record {
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
                "ID".to_string() => "ENSG00000235969.1".to_string(),
                "gene_id".to_string() => "ENSG00000235969.1".to_string(),
                "gene_type".to_string() => "unprocessed_pseudogene".to_string(),
                "gene_name".to_string() => "CHEK2P4".to_string(),
                "level".to_string() => "2".to_string(),
                "hgnc_id".to_string() => "HGNC:43580".to_string(),
                "havana_gene".to_string() => "OTTHUMG00000140400.1".to_string(),
            },
        };
        static ref CHEK2_EXON12: Gff3Record = Gff3Record {
            seq_id: "chr22".to_string(),
            source: "HAVANA".to_string(),
            feature_type: FeatureType::CDS,
            record_type: "CDS".to_string(),
            start: 28695710,
            end: 28695873,
            score: ".".to_string(),
            strand: GeneStrand::Reverse,
            phase: Some(0),
            attributes: indexmap! {
                "ID".to_string() => "CDS:ENST00000382580.6".to_string(),
                "Parent".to_string() => "ENST00000382580.6".to_string(),
                "gene_id".to_string() => "ENSG00000183765.22".to_string(),
                "transcript_id".to_string() => "ENST00000382580.6".to_string(),
                "gene_type".to_string() => "protein_coding".to_string(),
                "gene_name".to_string() => "CHEK2".to_string(),
                "transcript_type".to_string() => "protein_coding".to_string(),
                "transcript_name".to_string() => "CHEK2-203".to_string(),
                "exon_number".to_string() => "12".to_string(),
                "exon_id".to_string() => "ENSE00003559542.1".to_string(),
                "level".to_string() => "2".to_string(),
                "protein_id".to_string() => "ENSP00000372023.2".to_string(),
                "transcript_support_level".to_string() => "1".to_string(),
                "hgnc_id".to_string() => "HGNC:16627".to_string(),
                "tag".to_string() => "basic,appris_alternative_2,CCDS".to_string(),
                "ccdsid".to_string() => "CCDS33629.1".to_string(),
                "havana_gene".to_string() => "OTTHUMG00000151023.24".to_string(),
                "havana_transcript".to_string() => "OTTHUMT00000321016.2".to_string(),
            },
        };
    }

    const CHEK2P4_GENE_LINE: &str = "chr22\tHAVANA\tgene\t16503304\t16509475\t.\t+\t.\tID=ENSG00000235969.1;gene_id=ENSG00000235969.1;gene_type=unprocessed_pseudogene;gene_name=CHEK2P4;level=2;hgnc_id=HGNC:43580;havana_gene=OTTHUMG00000140400.1\n";
    const CHEK2_EXON12_LINE: &str = "chr22\tHAVANA\tCDS\t28695710\t28695873\t.\t-\t0\tID=CDS:ENST00000382580.6;Parent=ENST00000382580.6;gene_id=ENSG00000183765.22;transcript_id=ENST00000382580.6;gene_type=protein_coding;gene_name=CHEK2;transcript_type=protein_coding;transcript_name=CHEK2-203;exon_number=12;exon_id=ENSE00003559542.1;level=2;protein_id=ENSP00000372023.2;transcript_support_level=1;hgnc_id=HGNC:16627;tag=basic,appris_alternative_2,CCDS;ccdsid=CCDS33629.1;havana_gene=OTTHUMG00000151023.24;havana_transcript=OTTHUMT00000321016.2\n";

    #[allow(clippy::unreadable_literal)]
    #[test]
    fn test_gff3_parse() {
        assert_eq!(
            CHEK2P4_GENE_LINE.parse::<Gff3Record>().unwrap(),
            *CHEK2P4_GENE
        );

        assert_eq!(
            CHEK2_EXON12_LINE.parse::<Gff3Record>().unwrap(),
            *CHEK2_EXON12
        );
    }

    #[allow(clippy::unreadable_literal)]
    #[test]
    fn test_gff3_to_string() {
        assert_eq!(CHEK2P4_GENE_LINE, CHEK2P4_GENE.to_string());

        assert_eq!(CHEK2_EXON12_LINE, CHEK2_EXON12.to_string());
    }

    #[test]
    fn test_gff3_reader() {
        let test_data =
            autocompress::Decoder::suggest(
                &include_bytes!(
                    "../../testfiles/GENCODE/gencode.v33.basic.annotation.chr22.gff3.xz"
                )[..],
            )
            .unwrap();
        let test_reader = io::BufReader::new(test_data);
        let gff3_reader = Gff3Reader::new(test_reader);
        let values: Vec<_> = gff3_reader.map(|x| x.unwrap()).collect();
        assert_eq!(values[350], *CHEK2P4_GENE);
        assert_eq!(values[13123], *CHEK2_EXON12);
        assert_eq!(values.len(), 37856);
    }

    #[test]
    fn test_gff3_grouped_reader() {
        let test_data =
            autocompress::Decoder::suggest(
                &include_bytes!(
                    "../../testfiles/GENCODE/gencode.v33.basic.annotation.chr22.gff3.xz"
                )[..],
            )
            .unwrap();
        let test_reader = io::BufReader::new(test_data);
        let gff3_reader = Gff3Reader::new(test_reader);
        let grouped_gff3_reader = Gff3GroupedReader::new(gff3_reader);
        let values: Vec<_> = grouped_gff3_reader.map(|x| x.unwrap()).collect();

        assert_eq!(values[48].original_record, *CHEK2P4_GENE);
        assert_eq!(values.len(), 1388);
        // TODO: Check more details
    }
}
