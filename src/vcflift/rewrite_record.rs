use super::VCFHeaderRewriteTarget;
use crate::vcfparse::{
    CompleteVCFRecord, PartialVCFRecord, VCFParseError, VCFParseErrorKind, VCFRecord,
};
use log::warn;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::Write;
use std::str;

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::digit1;
use nom::combinator::opt;
use nom::multi::many1;
use nom::sequence::pair;

pub fn rewrite_info(
    record: &mut CompleteVCFRecord,
    rewrite_target: &VCFHeaderRewriteTarget,
    to_delete_index: &HashSet<usize>,
) {
    for (k, v) in record.info.iter_mut() {
        if rewrite_target.info_ref.contains(k.as_ref()) {
            // Reorder Number=R INFO
            let original_ref_data = v.remove(0);
            let new_ref_data = v
                .iter()
                .enumerate()
                .filter(|(i, _)| to_delete_index.contains(i))
                .map(|(_, x)| x.clone())
                .next()
                .unwrap_or_else(|| Cow::Owned(b".".to_vec()));
            *v = [new_ref_data]
                .iter()
                .cloned()
                .chain(
                    v.iter()
                        .enumerate()
                        .filter(|(i, _)| !to_delete_index.contains(i))
                        .map(|(_, x)| x.clone()),
                )
                .chain([original_ref_data].iter().cloned())
                .collect();
        } else if rewrite_target.info_alt.contains(k.as_ref()) {
            // Reorder Number=A INFO
            *v = v
                .iter()
                .enumerate()
                .filter(|(i, _)| !to_delete_index.contains(i))
                .map(|(_, x)| x.clone())
                .collect();
            v.push(Cow::Owned(b".".to_vec()));
        } else if rewrite_target.info_genotype.contains(k.as_ref()) {
            // TODO: implement here
        }
    }
}

pub fn rewrite_allele_frequency<'a>(
    original_record: &PartialVCFRecord<'a>,
    record: &mut CompleteVCFRecord,
    rewrite_target: &VCFHeaderRewriteTarget,
    to_delete_index: &HashSet<usize>,
) -> Result<(), VCFParseError> {
    for (k, v) in record.info.iter_mut() {
        if rewrite_target.allele_frequency.contains(k.as_ref()) {
            // Reorder Number=A INFO
            let mut frequency_sum = 0.0;
            for one in v.iter() {
                if one.as_ref() != &b"."[..] {
                    frequency_sum += str::from_utf8(one)?.parse::<f64>().map_err(|_| {
                        VCFParseErrorKind::FrequencyIsNotNumber {
                            line: original_record.line,
                        }
                    })?;
                }
            }

            let ref_freq = if frequency_sum > 1.0 {
                if frequency_sum > 1.01 {
                    warn!(
                        "sum of frequency is larger than 1: {} at {}:{} line {}",
                        frequency_sum,
                        str::from_utf8(original_record.contig()).unwrap(),
                        original_record.position(),
                        original_record.line
                    );
                }
                0.0
            } else {
                1.0 - frequency_sum
            };

            *v = v
                .iter()
                .enumerate()
                .filter(|(i, _)| !to_delete_index.contains(i))
                .map(|(_, x)| x.clone())
                .collect();
            v.push(Cow::Owned(format!("{}", ref_freq).as_bytes().to_vec()));
        }
    }

    Ok(())
}

pub fn rewrite_format(
    record: &mut CompleteVCFRecord,
    rewrite_target: &VCFHeaderRewriteTarget,
    to_delete_index: &HashSet<usize>,
) {
    for one_sample in record.call.iter_mut() {
        for (k, v) in record.format.iter().zip(one_sample.iter_mut()) {
            if rewrite_target.format_ref.contains(k.as_ref()) {
                // Reorder Number=R INFO
                let original_ref_data = v.remove(0);
                let new_ref_data = v
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| to_delete_index.contains(i))
                    .map(|(_, x)| x.clone())
                    .next()
                    .unwrap_or_else(|| Cow::Owned(b".".to_vec()));
                *v = [new_ref_data]
                    .iter()
                    .cloned()
                    .chain(
                        v.iter()
                            .enumerate()
                            .filter(|(i, _)| !to_delete_index.contains(i))
                            .map(|(_, x)| x.clone()),
                    )
                    .chain([original_ref_data].iter().cloned())
                    .collect();
            } else if rewrite_target.format_alt.contains(k.as_ref()) {
                // Reorder Number=A INFO
                *v = v
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| !to_delete_index.contains(i))
                    .map(|(_, x)| x.clone())
                    .collect();
                v.push(Cow::Owned(b".".to_vec()));
            } else {
                // TODO: implement here
            }
        }
    }
}

pub fn rewrite_gt<'a>(
    original_record: &PartialVCFRecord<'a>,
    record: &mut CompleteVCFRecord,
    _rewrite_target: &VCFHeaderRewriteTarget,
    to_delete_index: &HashSet<usize>,
) -> Result<(), VCFParseError> {
    let mut index_rewrite: HashMap<_, _> = (0..original_record.alternative().len())
        .filter(|x| !to_delete_index.contains(&x))
        .enumerate()
        .map(|(x, y)| (y + 1, x + 1))
        .collect();
    for i in to_delete_index.iter() {
        index_rewrite.insert(i + 1, 0);
    }
    index_rewrite.insert(0, record.alternative().len());
    //println!("Index rewrite: {:?}", index_rewrite);

    for one_sample in record.call.iter_mut() {
        for (k, v) in record.format.iter().zip(one_sample.iter_mut()) {
            if k == &&b"GT"[..] {
                if v.len() != 1 {
                    return Err(VCFParseErrorKind::InvalidGTRecord.into());
                }
                //println!("Parse GT: {:?}", str::from_utf8(&v[0]).unwrap());
                let parsed_gt: Vec<_> = parse_gt(&v[0])?
                    .iter()
                    .map(|x| (index_rewrite[&x.0], x.1))
                    .collect();
                *v = vec![Cow::Owned(format_gt(&parsed_gt).as_bytes().to_vec())];
            }
        }
    }

    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum GTSeparator {
    Phased,
    Unphased,
    None,
}

fn parse_gt(input: &[u8]) -> Result<Vec<(usize, GTSeparator)>, VCFParseError> {
    let (input, values) = many1(pair(digit1, opt(alt((tag(b"|"), tag(b"/"))))))(input)?;
    if !input.is_empty() {
        println!("Not empty input: {:?}", input);
        return Err(VCFParseErrorKind::InvalidGTRecord.into());
    }

    Ok(values
        .iter()
        .try_fold(Vec::new(), |mut acc, x| -> Result<_, VCFParseError> {
            let val: usize = str::from_utf8(x.0)?.parse()?;
            acc.push((
                val,
                match x.1 {
                    Some(b"/") => GTSeparator::Unphased,
                    Some(b"|") => GTSeparator::Phased,
                    None => GTSeparator::None,
                    _ => unreachable!(),
                },
            ));
            Ok(acc)
        })?)
}

fn format_gt(gt: &[(usize, GTSeparator)]) -> String {
    let mut result = String::new();
    for one in gt {
        write!(result, "{}", one.0).unwrap();
        match one.1 {
            GTSeparator::Phased => write!(result, "|").unwrap(),
            GTSeparator::Unphased => write!(result, "/").unwrap(),
            _ => (),
        }
    }
    result
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_gt() -> Result<(), VCFParseError> {
        assert_eq!(parse_gt(b"1")?, vec![(1, GTSeparator::None)]);
        assert_eq!(
            parse_gt(b"1/0")?,
            vec![(1, GTSeparator::Unphased), (0, GTSeparator::None)]
        );
        assert_eq!(
            parse_gt(b"1|0/2")?,
            vec![
                (1, GTSeparator::Phased),
                (0, GTSeparator::Unphased),
                (2, GTSeparator::None)
            ]
        );

        Ok(())
    }
}
