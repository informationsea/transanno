use super::VCFHeaderRewriteTarget;
use crate::vcfparse::{CompleteVCFRecord, PartialVCFRecord, VCFParseError, VCFRecord};
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
    let allele_number: HashMap<_, _> = record
        .info
        .iter()
        .filter(|(k, _)| rewrite_target.allele_number.contains(k.as_ref()))
        .cloned()
        .collect();
    for (k, v) in record.info.iter_mut() {
        if rewrite_target.allele_frequency.contains(k.as_ref()) {
            // Rewrite AF
            let frequency_sum = v
                .iter()
                .filter(|x| x.as_ref() != b".")
                .filter(|x| x.as_ref() != b"nan")
                .try_fold::<_, _, Result<f64, VCFParseError>>(0.0, |s, x| {
                    Ok(s + str::from_utf8(x)?.parse::<f64>().map_err(|_| {
                        VCFParseError::FrequencyIsNotNumber {
                            line: original_record.line,
                        }
                    })?)
                })?;

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
        } else if rewrite_target.allele_count.contains(k.as_ref()) {
            // Rewrite AC
            if let Some(corresponding_number) =
                rewrite_target.allele_count_to_allele_number.get(k.as_ref())
            {
                let corresponding_number: &[u8] = corresponding_number.as_ref();
                if let Some(number) = allele_number.get(corresponding_number) {
                    // println!(
                    //     "corresponding number: {}",
                    //     str::from_utf8(number[0].as_ref())?
                    // );
                    let number = str::from_utf8(number[0].as_ref())?.parse::<u64>()?;
                    let allele_count_sum = v
                        .iter()
                        .filter(|x| x.as_ref() != b".")
                        .try_fold::<_, _, Result<u64, VCFParseError>>(0, |s, x| {
                            Ok(s + str::from_utf8(x)?.parse::<u64>()?)
                        })?;
                    let ref_count = if number < allele_count_sum {
                        0
                    } else {
                        number - allele_count_sum
                    };

                    *v = v
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| !to_delete_index.contains(i))
                        .map(|(_, x)| x.clone())
                        .collect();
                    v.push(Cow::Owned(format!("{}", ref_count).as_bytes().to_vec()));
                }
            }
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
                    return Err(VCFParseError::InvalidGTRecord);
                }
                //println!("Parse GT: {:?}", str::from_utf8(&v[0]).unwrap());
                let mut parsed_gt: Vec<_> = parse_gt(&v[0])?
                    .iter()
                    .map(|x| (x.0.as_ref().map(|y| index_rewrite[y]), x.1))
                    .collect();
                if !parsed_gt.is_empty() && parsed_gt[0].1 == GTSeparator::Unphased {
                    parsed_gt.sort_by_key(|x| x.0);
                    let l = parsed_gt.len();
                    for (i, x) in parsed_gt.iter_mut().enumerate() {
                        if i + 1 == l {
                            x.1 = GTSeparator::None;
                        } else {
                            x.1 = GTSeparator::Unphased;
                        }
                    }
                }
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

fn parse_gt(input: &[u8]) -> Result<Vec<(Option<usize>, GTSeparator)>, VCFParseError> {
    let (input, values) = many1(pair(
        alt((digit1, tag(b"."))),
        opt(alt((tag(b"|"), tag(b"/")))),
    ))(input)?;
    if !input.is_empty() {
        println!("Not empty input: {:?}", input);
        return Err(VCFParseError::InvalidGTRecord);
    }

    Ok(values
        .iter()
        .try_fold(Vec::new(), |mut acc, x| -> Result<_, VCFParseError> {
            let val: Option<usize> = if x.0 == b"." {
                None
            } else {
                Some(str::from_utf8(x.0)?.parse()?)
            };
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

fn format_gt(gt: &[(Option<usize>, GTSeparator)]) -> String {
    let mut result = String::new();
    for one in gt {
        if let Some(x) = one.0 {
            write!(result, "{}", x).unwrap();
        } else {
            write!(result, ".").unwrap();
        }
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
        assert_eq!(parse_gt(b"1")?, vec![(Some(1), GTSeparator::None)]);
        assert_eq!(parse_gt(b".")?, vec![(None, GTSeparator::None)]);
        assert_eq!(
            parse_gt(b"1/0")?,
            vec![
                (Some(1), GTSeparator::Unphased),
                (Some(0), GTSeparator::None)
            ]
        );
        assert_eq!(
            parse_gt(b"1|./2")?,
            vec![
                (Some(1), GTSeparator::Phased),
                (None, GTSeparator::Unphased),
                (Some(2), GTSeparator::None)
            ]
        );

        Ok(())
    }
}
