mod rewrite_record;

use crate::defs::GenomeSequence;
use crate::variantlift::error::VariantLiftOverError;
use crate::variantlift::{LiftedVariant, VariantLiftOver};
use crate::vcfparse::{
    CompleteVCFRecord, PartialVCFRecord, VCFHeader, VCFHeaderItem, VCFParseError, VCFReader,
    VCFRecord, VCFWriter,
};
use crate::{chromosome_priority, LiftOverError, LiftOverErrorKind, Variant};
use log::{info, warn};
use regex::Regex;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::io::{self, prelude::*};
use std::str;

lazy_static! {
    pub static ref LIFT_SUCCESS_VCF_HEADER: Vec<VCFHeaderItem> = vec![
        VCFHeaderItem::parse(b"##INFO=<ID=MULTIMAP,Number=1,Type=Integer,Description=\"# of multi-mapped regions\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=REF_CHANGED,Number=0,Type=Flag,Description=\"Reference sequence is changed\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=ORIGINAL_REF,Number=1,Type=String,Description=\"Original reference sequence\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=ORIGINAL_CHROM,Number=1,Type=String,Description=\"Original chromosome\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=ORIGINAL_POS,Number=1,Type=Integer,Description=\"Original position\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=ORIGINAL_STRAND,Number=1,Type=String,Description=\"Original strand\">", 0).unwrap(),
    ];

    pub static ref LIFT_FAILED_VCF_HEADER: Vec<VCFHeaderItem> = vec![
        VCFHeaderItem::parse(b"##INFO=<ID=MULTIMAP,Number=1,Type=Integer,Description=\"# of multi-mapped regions\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=TESTED_CHROM,Number=1,Type=String,Description=\"Tested chromosome\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=TESTED_START,Number=1,Type=String,Description=\"Tested start\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=TESTED_END,Number=1,Type=String,Description=\"Tested end\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=FAILED_REASON,Number=1,Type=String,Description=\"Reason of liftOver failure\">", 0).unwrap(),
        VCFHeaderItem::parse(b"##INFO=<ID=PARTIAL_SUCCESS,Number=0,Type=Flag,Description=\"Variants in other tried region are succeeded to lift over\">", 0).unwrap(),
    ];

    pub static ref ALLELE_FREQUENCY_MATCH: Regex = Regex::new("^(.+_)?AF(_.+)?$").unwrap();
    pub static ref ALLELE_NUMBER_MATCH: Regex = Regex::new("^(.+_)?AN(_.+)?$").unwrap();
    pub static ref ALLELE_COUNT_MATCH: Regex = Regex::new("^(.+_)?AC(_.+)?$").unwrap();
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct VCFHeaderRewriteTarget {
    pub info_ref: HashSet<Vec<u8>>,
    pub info_alt: HashSet<Vec<u8>>,
    pub info_genotype: HashSet<Vec<u8>>,
    pub allele_frequency: HashSet<Vec<u8>>,
    pub allele_count: HashSet<Vec<u8>>,
    pub allele_count_to_allele_number: HashMap<Vec<u8>, Vec<u8>>,
    pub allele_number: HashSet<Vec<u8>>,
    pub format_ref: HashSet<Vec<u8>>,
    pub format_alt: HashSet<Vec<u8>>,
    pub format_genotype: HashSet<Vec<u8>>,
    pub format_gt: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum VCFRecordWrapper<'a> {
    Partial(PartialVCFRecord<'a>),
    Complete(CompleteVCFRecord<'a>),
}

impl<'a> VCFRecord for VCFRecordWrapper<'a> {
    fn contig(&self) -> &[u8] {
        match self {
            VCFRecordWrapper::Complete(c) => c.contig(),
            VCFRecordWrapper::Partial(p) => p.contig(),
        }
    }
    fn position(&self) -> u64 {
        match self {
            VCFRecordWrapper::Complete(c) => c.position(),
            VCFRecordWrapper::Partial(p) => p.position(),
        }
    }
    fn id(&self) -> &[u8] {
        match self {
            VCFRecordWrapper::Complete(c) => c.id(),
            VCFRecordWrapper::Partial(p) => p.id(),
        }
    }
    fn reference(&self) -> &[u8] {
        match self {
            VCFRecordWrapper::Complete(c) => c.reference(),
            VCFRecordWrapper::Partial(p) => p.reference(),
        }
    }
    fn alternative(&self) -> &[Cow<[u8]>] {
        match self {
            VCFRecordWrapper::Complete(c) => c.alternative(),
            VCFRecordWrapper::Partial(p) => p.alternative(),
        }
    }
    fn qual(&self) -> &[u8] {
        match self {
            VCFRecordWrapper::Complete(c) => c.qual(),
            VCFRecordWrapper::Partial(p) => p.qual(),
        }
    }
    fn filter(&self) -> &[u8] {
        match self {
            VCFRecordWrapper::Complete(c) => c.filter(),
            VCFRecordWrapper::Partial(p) => p.filter(),
        }
    }
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        match self {
            VCFRecordWrapper::Complete(c) => c.write(writer),
            VCFRecordWrapper::Partial(p) => p.write(writer),
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum VCFLiftOverResult<'a> {
    Succeeded(Vec<VCFRecordWrapper<'a>>),
    Failed(Box<PartialVCFRecord<'a>>),
}

#[derive(Debug, Copy, Clone, Hash, PartialEq, PartialOrd, Eq, Ord)]
pub struct VCFLiftOverParameters {
    pub allow_multimap: bool,
    pub acceptable_deletion: u64,
    pub acceptable_insertion: u64,
    pub do_not_rewrite_info: bool,
    pub do_not_rewrite_format: bool,
    pub do_not_rewrite_gt: bool,
    pub do_not_rewrite_allele_frequency: bool,
    pub do_not_rewrite_allele_count: bool,
    pub do_not_swap_ref_alt: bool,
    pub do_not_left_align_chain_file: bool,
    pub do_not_use_dot_when_alt_equal_to_ref: bool,
    pub do_not_prefer_cis_contig_when_multimap: bool,
}

impl VCFLiftOverParameters {
    pub fn new() -> Self {
        VCFLiftOverParameters {
            allow_multimap: false,
            acceptable_deletion: 3,
            acceptable_insertion: 3,
            do_not_rewrite_info: false,
            do_not_rewrite_format: false,
            do_not_rewrite_gt: false,
            do_not_rewrite_allele_frequency: false,
            do_not_rewrite_allele_count: false,
            do_not_swap_ref_alt: false,
            do_not_left_align_chain_file: false,
            do_not_use_dot_when_alt_equal_to_ref: false,
            do_not_prefer_cis_contig_when_multimap: false,
        }
    }

    pub fn allow_multimap(mut self, allow_multimap: bool) -> Self {
        self.allow_multimap = allow_multimap;
        self
    }

    pub fn acceptable_deletion(mut self, acceptable_deletion: u64) -> Self {
        self.acceptable_deletion = acceptable_deletion;
        self
    }

    pub fn acceptable_insertion(mut self, acceptable_insertion: u64) -> Self {
        self.acceptable_insertion = acceptable_insertion;
        self
    }

    pub fn do_not_rewrite_info(mut self, do_not_rewrite_info: bool) -> Self {
        self.do_not_rewrite_info = do_not_rewrite_info;
        self
    }

    pub fn do_not_rewrite_format(mut self, do_not_rewrite_format: bool) -> Self {
        self.do_not_rewrite_format = do_not_rewrite_format;
        self
    }

    pub fn do_not_rewrite_gt(mut self, do_not_rewrite_gt: bool) -> Self {
        self.do_not_rewrite_gt = do_not_rewrite_gt;
        self
    }

    pub fn do_not_rewrite_allele_frequency(mut self, do_not_rewrite_frequency: bool) -> Self {
        self.do_not_rewrite_allele_frequency = do_not_rewrite_frequency;
        self
    }

    pub fn do_not_rewrite_allele_count(mut self, do_not_rewrite_count: bool) -> Self {
        self.do_not_rewrite_allele_count = do_not_rewrite_count;
        self
    }

    pub fn do_not_swap_ref_alt(mut self, do_not_swap_ref_alt: bool) -> Self {
        self.do_not_swap_ref_alt = do_not_swap_ref_alt;
        self
    }

    pub fn do_not_left_align_chain_file(mut self, do_not_left_align_chain_file: bool) -> Self {
        self.do_not_left_align_chain_file = do_not_left_align_chain_file;
        self
    }

    pub fn do_not_use_dot_when_alt_equal_to_ref(
        mut self,
        do_not_use_dot_when_alt_equal_to_ref: bool,
    ) -> Self {
        self.do_not_use_dot_when_alt_equal_to_ref = do_not_use_dot_when_alt_equal_to_ref;
        self
    }

    pub fn do_not_prefer_cis_contig_when_multimap(
        mut self,
        do_not_prefer_cis_contig_when_multimap: bool,
    ) -> Self {
        self.do_not_prefer_cis_contig_when_multimap = do_not_prefer_cis_contig_when_multimap;
        self
    }
}

impl Default for VCFLiftOverParameters {
    fn default() -> Self {
        VCFLiftOverParameters::new()
    }
}

#[derive(Debug)]
pub struct VCFLiftOver<G: GenomeSequence> {
    variant_lift: VariantLiftOver<G>,
    param: VCFLiftOverParameters,
}

fn lift_header_show_warnings(
    rewrite_target: &mut VCFHeaderRewriteTarget,
    param: &VCFLiftOverParameters,
) {
    if !rewrite_target.info_genotype.is_empty()
        && !param.do_not_swap_ref_alt
        && !param.do_not_rewrite_info
    {
        warn!("Rewriting INFO with Number=G is not implemented.");
    }

    if !rewrite_target.format_genotype.is_empty()
        && !param.do_not_swap_ref_alt
        && !param.do_not_rewrite_format
    {
        warn!("Rewriting FORMAT with Number=G is not implemented.");
    }
}

impl<G: GenomeSequence> VCFLiftOver<G> {
    pub fn new(variant_lift: VariantLiftOver<G>, param: VCFLiftOverParameters) -> Self {
        VCFLiftOver {
            variant_lift,
            param,
        }
    }

    pub fn lift_header(
        &self,
        header: &VCFHeader,
    ) -> Result<(VCFHeader, VCFHeaderRewriteTarget), LiftOverError> {
        let mut new_header_items = Vec::new();
        let mut rewrite_target = VCFHeaderRewriteTarget::default();
        let mut allele_count: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        let mut allele_number: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();

        for one_item in &header.header_items {
            if one_item.key == b"contig" {
                self.check_contig_length_helper(&one_item)?;
            } else {
                new_header_items.push(one_item.clone())
            }

            if one_item.key == b"source" {
                if one_item.value == b"ClinVar" && !self.param.do_not_swap_ref_alt {
                    warn!("Input VCF looks like ClinVar. You should add --noswap option to avoid swapping reference and alternative allele");
                }
                if one_item.value.starts_with(b"COSMIC") && !self.param.do_not_swap_ref_alt {
                    warn!("Input VCF looks like COSMIC. You should add --noswap option to avoid swapping reference and alternative allele");
                }
            }

            if !self.param.do_not_rewrite_info && one_item.key == b"INFO" {
                if let Some(id) = one_item.detail.get(&b"ID"[..]) {
                    if let Some(number) = one_item.detail.get(&b"Number"[..]) {
                        if !self.param.do_not_rewrite_allele_frequency
                            && ALLELE_FREQUENCY_MATCH.is_match(str::from_utf8(id).unwrap())
                            && number == b"A"
                        {
                            rewrite_target.allele_frequency.insert(id.to_vec());
                        } else if !self.param.do_not_rewrite_allele_count
                            && ALLELE_COUNT_MATCH.is_match(str::from_utf8(id).unwrap())
                            && number == b"A"
                        {
                            rewrite_target.allele_count.insert(id.to_vec());

                            let m = ALLELE_COUNT_MATCH
                                .captures(str::from_utf8(id).unwrap())
                                .unwrap();
                            let mut key: Vec<u8> = Vec::new();
                            key.extend_from_slice(
                                m.get(1).map(|x| x.as_str().as_bytes()).unwrap_or(b""),
                            );
                            key.push(b' ');
                            key.extend_from_slice(
                                m.get(2).map(|x| x.as_str().as_bytes()).unwrap_or(b""),
                            );
                            allele_count.insert(id.to_vec(), key);
                        } else if !self.param.do_not_rewrite_allele_count
                            && ALLELE_NUMBER_MATCH.is_match(str::from_utf8(id).unwrap())
                            && number == b"1"
                        {
                            rewrite_target.allele_number.insert(id.to_vec());
                            let m = ALLELE_NUMBER_MATCH
                                .captures(str::from_utf8(id).unwrap())
                                .unwrap();
                            let mut key: Vec<u8> = Vec::new();
                            key.extend_from_slice(
                                m.get(1).map(|x| x.as_str().as_bytes()).unwrap_or(b""),
                            );
                            key.push(b' ');
                            key.extend_from_slice(
                                m.get(2).map(|x| x.as_str().as_bytes()).unwrap_or(b""),
                            );
                            allele_number.insert(key, id.to_vec());
                        } else {
                            match &number as &[u8] {
                                b"G" => rewrite_target.info_genotype.insert(id.to_vec()),
                                b"R" => rewrite_target.info_ref.insert(id.to_vec()),
                                b"A" => rewrite_target.info_alt.insert(id.to_vec()),
                                _ => false,
                            };
                        }
                    }
                }
            }

            if !self.param.do_not_rewrite_format && one_item.key == b"FORMAT" {
                if let Some(id) = one_item.detail.get(&b"ID"[..]) {
                    if let Some(number) = one_item.detail.get(&b"Number"[..]) {
                        match &number as &[u8] {
                            b"G" => rewrite_target.format_genotype.insert(id.to_vec()),
                            b"R" => rewrite_target.format_ref.insert(id.to_vec()),
                            b"A" => rewrite_target.format_alt.insert(id.to_vec()),
                            _ => false,
                        };
                    }

                    if id == b"GT" {
                        rewrite_target.format_gt = true;
                    }
                }
            }
        }

        //println!("allele count: {:?}", allele_count);
        //println!("allele number: {:?}", allele_number);
        for (k, v) in allele_count.iter() {
            if let Some(x) = allele_number.get(v) {
                rewrite_target
                    .allele_count_to_allele_number
                    .insert(k.to_vec(), x.to_vec());
            }
        }

        lift_header_show_warnings(&mut rewrite_target, &self.param);

        let mut chromosomes = self
            .variant_lift
            .position_liftover()
            .query_chromosomes()
            .to_vec();
        chromosomes.sort_by_key(|x| (chromosome_priority(&x.name), x.length));
        for one_chromosome in chromosomes {
            new_header_items.push(VCFHeaderItem::parse(
                &format!(
                    "##contig=<ID={},length={}>",
                    one_chromosome.name, one_chromosome.length
                )
                .as_bytes(),
                0,
            )?);
        }

        for one_item in LIFT_SUCCESS_VCF_HEADER.iter() {
            new_header_items.push(one_item.clone());
        }
        Ok((
            VCFHeader {
                header_items: new_header_items,
                samples: header.samples.clone(),
            },
            rewrite_target,
        ))
    }

    pub fn lift_record<'a>(
        &mut self,
        record: &'a PartialVCFRecord,
        rewrite_target: &VCFHeaderRewriteTarget,
    ) -> Result<VCFLiftOverResult<'a>, LiftOverError> {
        let mut original_variant: Variant = record.into();

        // add chr prefix if it is required
        if self
            .variant_lift
            .position_liftover()
            .reference_chromosome_by_name(str::from_utf8(record.contig()).unwrap())
            .is_none()
            && !original_variant.chromosome.starts_with("chr")
        {
            let chr_added_name = format!("chr{}", original_variant.chromosome);
            if self
                .variant_lift
                .position_liftover()
                .reference_chromosome_by_name(&chr_added_name)
                .is_some()
            {
                original_variant.chromosome = chr_added_name;
            }
        }

        let lifted_variant = self.variant_lift.lift_variant(
            &original_variant,
            self.param.acceptable_deletion,
            self.param.acceptable_insertion,
        )?;

        let mut succeeded_records = Vec::new();

        let mut failed_reasons = HashSet::new();
        for one in lifted_variant.iter() {
            match one {
                Ok(ok) => {
                    let new_record = merge_to_vcf(&ok, record, &self.param, rewrite_target)?;
                    succeeded_records.push(new_record);
                }
                Err(e) => match e {
                    VariantLiftOverError::UnacceptableLargeDeletion { .. } => {
                        failed_reasons.insert("UNACCEPTABLE_LARGE_DELETION");
                    }
                    VariantLiftOverError::UnacceptableLargeInsertion { .. } => {
                        failed_reasons.insert("UNACCEPTABLE_LARGE_INSERTION");
                    }
                    VariantLiftOverError::UnknownSequenceName(_chrom) => {
                        failed_reasons.insert("UNKNOWN_SEQUENCE_NAME");
                    }
                    VariantLiftOverError::ReferenceSequenceIsNotMatch => {
                        failed_reasons.insert("UNEXPECTED_REF");
                    }
                },
            }
        }

        // create failed record if succeeded record list is empty and some failed reasons are found.
        if succeeded_records.is_empty() && !failed_reasons.is_empty() {
            let mut new_record = record.clone();
            if new_record.unparsed_info == &b"."[..] {
                new_record.unparsed_info = Cow::Owned(b"FAILED_REASON=".to_vec());
            } else {
                write!(new_record.unparsed_info.to_mut(), ";FAILED_REASON=",)?;
            }
            for (i, one_reason) in failed_reasons.iter().enumerate() {
                if i != 0 {
                    write!(new_record.unparsed_info.to_mut(), ",")?;
                }
                write!(new_record.unparsed_info.to_mut(), "{}", one_reason)?;
            }
            return Ok(VCFLiftOverResult::Failed(Box::new(new_record)));
        }

        // use cis mapped variants if it was found.
        if !self.param.do_not_prefer_cis_contig_when_multimap
            && succeeded_records
                .iter()
                .any(|x| x.contig() != original_variant.chromosome.as_bytes())
            && succeeded_records
                .iter()
                .any(|x| x.contig() == original_variant.chromosome.as_bytes())
        {
            succeeded_records.retain(|x| x.contig() == original_variant.chromosome.as_bytes());
        }

        // multi-map check
        if succeeded_records.len() > 1 {
            if !self.param.allow_multimap {
                let mut new_record = record.clone();
                if new_record.unparsed_info == &b"."[..] {
                    new_record.unparsed_info = Cow::Borrowed(b"");
                } else {
                    new_record.unparsed_info.to_mut().push(b';');
                }
                write!(
                    new_record.unparsed_info.to_mut(),
                    "FAILED_REASON=MULTIMAP;MULTIMAP={}",
                    succeeded_records.len()
                )?;
                return Ok(VCFLiftOverResult::Failed(Box::new(new_record)));
            } else {
                let success_count = succeeded_records.len();
                succeeded_records
                    .iter_mut()
                    .for_each(|mut x| add_multimap_info_helper(&mut x, success_count));
            }
        }

        // no chain was found.
        if succeeded_records.is_empty() {
            let mut new_record = record.clone();
            write!(new_record.unparsed_info.to_mut(), "FAILED_REASON=NO_CHAIN")?;
            return Ok(VCFLiftOverResult::Failed(Box::new(new_record)));
        }

        Ok(VCFLiftOverResult::Succeeded(succeeded_records))
    }

    pub fn lift_vcf<R: Read, W1: Write, W2: Write>(
        &mut self,
        reader: R,
        success_writer: W1,
        failed_writer: W2,
    ) -> Result<(), LiftOverError> {
        let mut vcf_reader = VCFReader::new(reader)?;
        let lifted_header = self.lift_header(&vcf_reader.header)?;
        let mut success_vcf_writer =
            VCFWriter::new(io::BufWriter::new(success_writer), lifted_header.0)?;

        let mut failed_header = vcf_reader.header.clone();
        for one_item in LIFT_FAILED_VCF_HEADER.iter() {
            failed_header.header_items.push(one_item.clone());
        }

        let mut failed_vcf_writer = VCFWriter::new(failed_writer, failed_header)?;

        let mut multi_allelic_warned = (self.param.do_not_rewrite_allele_frequency
            || lifted_header.1.allele_frequency.is_empty())
            && (self.param.do_not_rewrite_allele_count || lifted_header.1.allele_count.is_empty());
        let mut last_position = 0;

        let mut succeeded_records = 0;
        let mut failed_records = 0;
        while let Some(original_record) = vcf_reader.next_record()? {
            if !multi_allelic_warned && last_position == original_record.position {
                warn!("Multi allelic sites should be merged into one line to rewrite allele frequency/count correctly.");
                warn!("Please merge multi allelic sites with `bcftools norm -m +any` command");
                multi_allelic_warned = true;
            }
            last_position = original_record.position;

            let lifted_record = self.lift_record(&original_record, &lifted_header.1)?;

            match lifted_record {
                VCFLiftOverResult::Succeeded(succeeded) => {
                    for one_success in succeeded {
                        success_vcf_writer.write_record(&one_success)?;
                    }
                    succeeded_records += 1;
                }
                VCFLiftOverResult::Failed(failed) => {
                    failed_vcf_writer.write_record(failed.as_ref())?;
                    failed_records += 1;
                }
            }

            if (succeeded_records + failed_records) % 1_000_000 == 0 {
                info!(
                    "Processed {} entries at {}:{}",
                    succeeded_records + failed_records,
                    str::from_utf8(&original_record.contig).unwrap(),
                    original_record.position
                );
            }
        }

        eprintln!("    Total record: {}", succeeded_records + failed_records);
        eprintln!("   Mapped record: {}", succeeded_records);
        eprintln!(" Unmapped record: {}", failed_records);

        Ok(())
    }

    fn check_contig_length_helper(&self, one_item: &VCFHeaderItem) -> Result<(), LiftOverError> {
        if let Some(length) = one_item.detail.get(&b"length"[..]) {
            if let Some(id) = one_item.detail.get(&b"ID"[..]) {
                if let Ok(length) = str::from_utf8(length).unwrap().parse::<u64>() {
                    if let Some(known_chromosome) = self
                        .variant_lift
                        .position_liftover()
                        .reference_chromosome_by_name(str::from_utf8(id).unwrap())
                    {
                        if length != known_chromosome.length {
                            return Err(LiftOverErrorKind::ReferenceChromosomeLengthIsNotMatch(
                                str::from_utf8(id).unwrap().to_string(),
                            )
                            .into());
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

fn add_multimap_info_helper(record: &mut VCFRecordWrapper, count: usize) {
    match record {
        VCFRecordWrapper::Partial(ref mut p) => {
            if p.unparsed_info == &b"."[..] {
                p.unparsed_info = Cow::Borrowed(b"");
            } else {
                p.unparsed_info.to_mut().push(b';');
            }
            write!(p.unparsed_info.to_mut(), "MULTIMAP={}", count).unwrap();
        }
        VCFRecordWrapper::Complete(ref mut c) => {
            c.info.push((
                Cow::Owned(b"MULTIMAP".to_vec()),
                vec![Cow::Owned(format!("{}", count).as_bytes().to_vec())],
            ));
        }
    }
}

fn merge_to_vcf<'a>(
    variant: &LiftedVariant,
    record: &PartialVCFRecord<'a>,
    param: &VCFLiftOverParameters,
    rewrite_target: &VCFHeaderRewriteTarget,
) -> Result<VCFRecordWrapper<'a>, VCFParseError> {
    if variant.reference_changed && !param.do_not_swap_ref_alt {
        Ok(VCFRecordWrapper::Complete(swap_ref_alt(
            variant,
            record,
            param,
            rewrite_target,
        )?))
    } else {
        Ok(VCFRecordWrapper::Partial(merge_to_vcf_simple(
            variant, record, param,
        )))
    }
}

fn merge_to_vcf_simple<'a>(
    variant: &LiftedVariant,
    record: &PartialVCFRecord<'a>,
    param: &VCFLiftOverParameters,
) -> PartialVCFRecord<'a> {
    let mut new_record = record.clone();

    // update contig/ref/alt
    new_record.contig = Cow::Owned(variant.chromosome.as_bytes().to_vec());
    new_record.reference = Cow::Owned(variant.reference.to_vec());
    new_record.alternative = variant
        .alternative
        .iter()
        .map(|x| {
            if !param.do_not_use_dot_when_alt_equal_to_ref && x == &variant.reference {
                Cow::Owned(b".".to_vec())
            } else {
                Cow::Owned(x.to_vec())
            }
        })
        .collect();
    new_record.position = variant.position + 1;

    let new_info = new_record.unparsed_info.to_mut();
    if new_info as &[u8] == &b"."[..] {
        new_info.clear();
    } else {
        new_info.push(b';');
    }
    new_info.append(&mut b"ORIGINAL_CHROM=".to_vec());
    new_info.append(&mut record.contig().to_vec());
    write!(
        new_info,
        ";ORIGINAL_POS={};ORIGINAL_STRAND={}",
        record.position(),
        variant.strand
    )
    .unwrap();

    if variant.reference_changed {
        new_info.append(&mut b";ORIGINAL_REF=".to_vec());
        new_info.append(&mut variant.original_reference.to_vec());
    }

    new_record
}

fn swap_ref_alt<'a>(
    lifted_variant: &LiftedVariant,
    original_record: &PartialVCFRecord<'a>,
    param: &VCFLiftOverParameters,
    rewrite_target: &VCFHeaderRewriteTarget,
) -> Result<CompleteVCFRecord<'a>, VCFParseError> {
    let mut record = original_record.clone().complete_parse()?;
    record.contig = Cow::Owned(lifted_variant.chromosome.as_bytes().to_vec());
    record.reference = Cow::Owned(lifted_variant.reference.to_vec());
    record.alternative = lifted_variant
        .alternative
        .iter()
        .map(|x| Cow::Owned(x.to_vec()))
        .collect();
    record.position = lifted_variant.position + 1;

    record.info.push((
        Cow::Owned(b"ORIGINAL_CHROM".to_vec()),
        vec![Cow::Owned(original_record.contig().to_vec())],
    ));
    record.info.push((
        Cow::Owned(b"ORIGINAL_POS".to_vec()),
        vec![Cow::Owned(
            format!("{}", original_record.position())
                .as_bytes()
                .to_vec(),
        )],
    ));
    record.info.push((
        Cow::Owned(b"ORIGINAL_STRAND".to_vec()),
        vec![Cow::Owned(
            format!("{}", lifted_variant.strand).as_bytes().to_vec(),
        )],
    ));

    if record.reference() == &lifted_variant.original_reference[..] {
        return Ok(record);
    }

    record
        .info
        .push((Cow::Owned(b"REF_CHANGED".to_vec()), vec![]));
    record.info.push((
        Cow::Owned(b"ORIGINAL_REF".to_vec()),
        vec![Cow::Owned(lifted_variant.original_reference.to_vec())],
    ));

    let to_delete_index: HashSet<_> = record
        .alternative()
        .iter()
        .enumerate()
        .filter(|(_, x)| x.as_ref() == record.reference())
        .map(|(i, _)| i)
        .collect();
    record.alternative = record
        .alternative
        .iter()
        .enumerate()
        .filter(|(i, _)| !to_delete_index.contains(i))
        .map(|(_, x)| x.clone())
        .chain(vec![Cow::Owned(lifted_variant.original_reference.to_vec())].into_iter())
        .collect();

    if !param.do_not_rewrite_info {
        rewrite_record::rewrite_info(&mut record, rewrite_target, &to_delete_index);
    }

    if !param.do_not_rewrite_allele_frequency {
        rewrite_record::rewrite_allele_frequency(
            original_record,
            &mut record,
            rewrite_target,
            &to_delete_index,
        )?;
    }

    if !param.do_not_rewrite_format {
        rewrite_record::rewrite_format(&mut record, rewrite_target, &to_delete_index);
    }

    if !param.do_not_rewrite_gt {
        rewrite_record::rewrite_gt(
            original_record,
            &mut record,
            rewrite_target,
            &to_delete_index,
        )?;
    }

    // TODO: add original chrom/pos/ref record

    Ok(record)
}

#[cfg(test)]
mod test;
