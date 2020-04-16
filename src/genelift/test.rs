use super::*;
use crate::poslift::PositionLiftOver;
use autocompress::{create, open};
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::fs;

#[derive(Debug, Clone, PartialEq)]
struct TestFeature {
    pub seq_id: String,
    pub feature_type: FeatureType,
    pub start: u64,
    pub end: u64,
    pub strand: GeneStrand,
    pub attributes: HashMap<String, String>,
}

impl Display for TestFeature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{:?}", self)
    }
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

fn load_hg38_to_hg19_lift() -> GeneLiftOver<PositionLiftOver> {
    let hg38_to_hg19_chain = open("testfiles/hg38ToHg19/hg38ToHg19.over.chain.chr22").unwrap();
    let region_lift = PositionLiftOver::load(hg38_to_hg19_chain).unwrap();
    GeneLiftOver::new(region_lift)
}

#[allow(clippy::unreadable_literal)]
#[test]
fn test_simple_feature() {
    let gene_lift = load_hg38_to_hg19_lift();

    // test exon
    let original = TestFeature {
        feature_type: FeatureType::Exon,
        seq_id: "chr22".to_string(),
        start: 21867322,
        end: 21867629,
        strand: GeneStrand::Forward,
        attributes: HashMap::new(),
    };

    let lifted_feature = gene_lift.lift_single_feature(&original).unwrap();
    assert_eq!(
        lifted_feature,
        vec![LiftedFeature {
            seq_id: "chr22".to_string(),
            start: 22221612,
            end: 22221919,
            strand: GeneStrand::Forward,
            feature: &original,
            changes: vec![RegionChangeOp::Aligned(22221919 - 22221612)],
        }]
    );

    assert_eq!(
        lifted_feature[0].apply(),
        TestFeature {
            feature_type: FeatureType::Exon,
            seq_id: "chr22".to_string(),
            start: 22221612,
            end: 22221919,
            strand: GeneStrand::Forward,
            attributes: vec![
                ("ORIGINAL_CHROM".to_string(), "chr22".to_string()),
                ("ORIGINAL_START".to_string(), "21867322".to_string()),
                ("ORIGINAL_END".to_string(), "21867629".to_string()),
                ("ORIGINAL_STRAND".to_string(), "+".to_string()),
                ("CIGER".to_string(), "M307".to_string())
            ]
            .into_iter()
            .collect(),
        }
    );

    // test transcript
    let original = Transcript {
        children: vec![TestFeature {
            feature_type: FeatureType::Exon,
            seq_id: "chr22".to_string(),
            start: 22822776,
            end: 22822871,
            strand: GeneStrand::Forward,
            attributes: HashMap::new(),
        }],
        original_record: TestFeature {
            feature_type: FeatureType::Transcript,
            seq_id: "chr22".to_string(),
            start: 22822776,
            end: 22822871,
            strand: GeneStrand::Forward,
            attributes: HashMap::new(),
        },
    };

    let lifted_transcript = gene_lift.lift_transcript_feature(&original).unwrap();

    assert_eq!(
        lifted_transcript,
        LiftedTranscript {
            features: vec![LiftedFeature {
                seq_id: "chr22".to_string(),
                start: 23165270,
                end: 23165365,
                changes: vec![RegionChangeOp::Aligned(23165365 - 23165270)],
                feature: &original.children[0],
                strand: GeneStrand::Forward,
            }],
            transcript: &original,
            seq_id: "chr22".to_string(),
            start: 23165270,
            end: 23165365,
            strand: GeneStrand::Forward,
        }
    );

    assert_eq!(
        lifted_transcript.apply(),
        Transcript {
            children: vec![TestFeature {
                feature_type: FeatureType::Exon,
                seq_id: "chr22".to_string(),
                start: 23165270,
                end: 23165365,
                strand: GeneStrand::Forward,
                attributes: vec![
                    ("ORIGINAL_CHROM".to_string(), "chr22".to_string()),
                    ("ORIGINAL_START".to_string(), "22822776".to_string()),
                    ("ORIGINAL_END".to_string(), "22822871".to_string()),
                    ("ORIGINAL_STRAND".to_string(), "+".to_string()),
                    ("CIGER".to_string(), "M95".to_string())
                ]
                .into_iter()
                .collect(),
            }],
            original_record: TestFeature {
                feature_type: FeatureType::Transcript,
                seq_id: "chr22".to_string(),
                start: 23165270,
                end: 23165365,
                strand: GeneStrand::Forward,
                attributes: vec![
                    ("ORIGINAL_CHROM".to_string(), "chr22".to_string()),
                    ("ORIGINAL_START".to_string(), "22822776".to_string()),
                    ("ORIGINAL_END".to_string(), "22822871".to_string()),
                    ("ORIGINAL_STRAND".to_string(), "+".to_string()),
                ]
                .into_iter()
                .collect(),
            },
        }
    );

    // test gene
    let original = Gene {
        transcripts: vec![Transcript {
            children: vec![TestFeature {
                feature_type: FeatureType::Exon,
                seq_id: "chr22".to_string(),
                start: 22822776,
                end: 22822871,
                strand: GeneStrand::Forward,
                attributes: HashMap::new(),
            }],
            original_record: TestFeature {
                feature_type: FeatureType::Transcript,
                seq_id: "chr22".to_string(),
                start: 22822776,
                end: 22822871,
                strand: GeneStrand::Forward,
                attributes: HashMap::new(),
            },
        }],
        original_record: TestFeature {
            feature_type: FeatureType::Gene,
            seq_id: "chr22".to_string(),
            start: 22822776,
            end: 22822871,
            strand: GeneStrand::Forward,
            attributes: HashMap::new(),
        },
    };

    assert_eq!(
        gene_lift.lift_gene_feature(&original).unwrap(),
        LiftedGene {
            transcripts: vec![LiftedTranscript {
                transcript: &original.transcripts[0],
                features: vec![LiftedFeature {
                    seq_id: "chr22".to_string(),
                    start: 23165270,
                    end: 23165365,
                    strand: GeneStrand::Forward,
                    changes: vec![RegionChangeOp::Aligned(23165365 - 23165270)],
                    feature: &original.transcripts[0].children[0],
                }],
                seq_id: "chr22".to_string(),
                start: 23165270,
                end: 23165365,
                strand: GeneStrand::Forward,
            }],
            failed_transcripts: vec![],
            seq_id: "chr22".to_string(),
            start: 23165270,
            end: 23165365,
            gene: &original,
            strand: GeneStrand::Forward
        }
    );
}

// TODO: test complex features

use crate::geneparse::gff3::{Gff3GroupedReader, Gff3Reader};
use crate::geneparse::gtf::{GtfGroupedReader, GtfReader};
use std::io::{BufReader, Write};

#[test]
fn test_gff3_real_check2() -> Result<(), FeatureLiftError> {
    let gff3_data =
        &include_bytes!("../../testfiles/gene/gencode.v30.basic.annotation.CHEK2-MCHR1.gff3")[..];
    let gff3_genes = Gff3GroupedReader::new(Gff3Reader::new(BufReader::new(gff3_data)));
    let gene_lift = load_hg38_to_hg19_lift();

    fs::create_dir_all("target/test-output/gene").unwrap();

    let mut writer =
        create("./target/test-output/gene/gencode.v30.basic.annotation.CHEK2-MCHR1-lifted.gff3")
            .unwrap();
    for one in gff3_genes {
        let one = one.unwrap();
        match gene_lift.lift_gene_feature(&one) {
            Ok(lifted_gene) => {
                let to_write = lifted_gene.apply();
                write!(writer, "{}", to_write).unwrap();
            }
            Err(e) => {
                println!(
                    "gene lift error: {:?} {:?}",
                    e.error,
                    e.gene.attribute("ID")
                );
                for one_failed_transcript in e.failed_transcripts {
                    println!(
                        "   transcript lift error: {:?} {:?} {} {}",
                        one_failed_transcript.error,
                        one_failed_transcript.transcript.attribute("ID"),
                        one_failed_transcript.failed_features.len(),
                        one_failed_transcript.features.len(),
                    );
                }
            }
        }
    }

    Ok(())
}

#[test]
fn test_gtf_real_check2() -> Result<(), FeatureLiftError> {
    fs::create_dir_all("target/test-output/gene").unwrap();

    let gtf_data =
        &include_bytes!("../../testfiles/gene/gencode.v31.annotation.CHEK2-MCHR1.gtf")[..];
    let gtf_genes = GtfGroupedReader::new(GtfReader::new(BufReader::new(gtf_data)));
    let gene_lift = load_hg38_to_hg19_lift();
    let mut writer =
        create("./target/test-output/gene/gencode.v31.annotation.CHEK2-MCHR1-lifted.gff3").unwrap();
    for one in gtf_genes {
        let one = one.unwrap();
        match gene_lift.lift_gene_feature(&one) {
            Ok(lifted_gene) => {
                let to_write = lifted_gene.apply();
                write!(writer, "{}", to_write).unwrap();
            }
            Err(e) => {
                println!(
                    "gene lift error: {:?} {:?}",
                    e.error,
                    e.gene.attribute("ID")
                );
                for one_failed_transcript in e.failed_transcripts {
                    println!(
                        "   transcript lift error: {:?} {:?} {} {}",
                        one_failed_transcript.error,
                        one_failed_transcript.transcript.attribute("ID"),
                        one_failed_transcript.failed_features.len(),
                        one_failed_transcript.features.len(),
                    );
                }
            }
        }
    }

    Ok(())
}
