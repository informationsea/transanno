use super::{as_record_error, CompleteVCFRecord, InfoPair, VCFParseError, VCFRecord};
use nom::bytes::complete::{tag, take_till, take_while};
use nom::character::complete::digit1;
use nom::multi::separated_list1;
use nom::IResult;
use std::borrow::Cow;
use std::fmt::{self, Debug};
use std::io::{self, Write};
use std::str;

#[derive(PartialEq, Clone)]
pub struct PartialVCFRecord<'a> {
    pub raw: &'a [u8],
    pub line: u32,
    pub contig: Cow<'a, [u8]>,
    pub position: u64,
    pub id: Cow<'a, [u8]>,
    pub reference: Cow<'a, [u8]>,
    pub alternative: Vec<Cow<'a, [u8]>>,
    pub qual: Cow<'a, [u8]>,
    pub filter: Cow<'a, [u8]>,
    pub unparsed_info: Cow<'a, [u8]>,
    pub(crate) original_unparsed_info: &'a [u8],
    pub other: &'a [u8],
}

impl<'a> Debug for PartialVCFRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(f, "PartialVCFRecord {{ ")?;
        write!(f, "raw: b\"{}\", ", str::from_utf8(self.raw).unwrap())?;
        write!(
            f,
            "contig: b\"{}\", ",
            str::from_utf8(&self.contig).unwrap()
        )?;
        write!(f, "position: b\"{}\", ", self.position)?;
        write!(f, "id: b\"{}\", ", str::from_utf8(&self.id).unwrap())?;
        write!(
            f,
            "reference: b\"{}\", ",
            str::from_utf8(&self.reference).unwrap()
        )?;
        write!(f, "alternative: vec![")?;
        for (i, x) in self.alternative.iter().enumerate() {
            if i != 0 {
                write!(f, ", ")?;
            }
            write!(f, "b\"{}\"", str::from_utf8(x).unwrap())?;
        }
        write!(f, "], ")?;
        write!(f, "qual: b\"{}\", ", str::from_utf8(&self.qual).unwrap())?;
        write!(
            f,
            "filter: b\"{}\", ",
            str::from_utf8(&self.filter).unwrap()
        )?;
        write!(
            f,
            "unparsed_info: b\"{}\", ",
            str::from_utf8(&self.unparsed_info).unwrap()
        )?;
        write!(f, "other: b\"{}\"}}", str::from_utf8(self.other).unwrap())?;

        Ok(())
    }
}

impl<'a> VCFRecord for PartialVCFRecord<'a> {
    fn contig(&self) -> &[u8] {
        &self.contig
    }
    fn position(&self) -> u64 {
        self.position
    }
    fn id(&self) -> &[u8] {
        &self.id
    }
    fn reference(&self) -> &[u8] {
        &self.reference
    }
    fn alternative(&self) -> &[Cow<[u8]>] {
        &self.alternative
    }
    fn qual(&self) -> &[u8] {
        &self.qual
    }
    fn filter(&self) -> &[u8] {
        &self.filter
    }
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.contig)?;
        writer.write_all(b"\t")?;
        write!(writer, "{}", self.position)?;
        writer.write_all(b"\t")?;
        writer.write_all(&self.id)?;
        writer.write_all(b"\t")?;
        writer.write_all(&self.reference)?;
        writer.write_all(b"\t")?;
        for (i, one) in self.alternative.iter().enumerate() {
            if i != 0 {
                writer.write_all(b",")?;
            }
            writer.write_all(one)?;
        }
        writer.write_all(b"\t")?;
        writer.write_all(if self.qual.is_empty() {
            b"."
        } else {
            &self.qual
        })?;
        writer.write_all(b"\t")?;
        writer.write_all(if self.filter.is_empty() {
            b"."
        } else {
            &self.filter
        })?;
        writer.write_all(b"\t")?;
        writer.write_all(if self.unparsed_info.is_empty() {
            b"."
        } else {
            &self.unparsed_info
        })?;
        if !self.other.is_empty() {
            writer.write_all(b"\t")?;
            writer.write_all(&self.other)?;
        }
        writer.write_all(b"\n")?;
        Ok(())
    }
}

impl<'a> PartialVCFRecord<'a> {
    pub fn parse_vcf(line: u32, input: &[u8]) -> Result<PartialVCFRecord, VCFParseError> {
        match parse_vcf_internal(line, input) {
            Ok(x) => Ok(x),
            Err(nom::Err::Failure(e)) | Err(nom::Err::Error(e)) => {
                Err(as_record_error(line, input, e))
            }
            _ => unreachable!(),
        }
    }

    pub fn complete_parse(self) -> Result<CompleteVCFRecord<'a>, VCFParseError> {
        let raw = self.raw;
        let line = self.line;
        match self.complete_parse_helper() {
            Ok(x) => Ok(x),
            Err(nom::Err::Failure(e)) | Err(nom::Err::Error(e)) => {
                Err(as_record_error(line, &raw, e))
            }
            _ => unreachable!(),
        }
    }

    fn complete_parse_helper(
        self,
    ) -> Result<CompleteVCFRecord<'a>, nom::Err<impl nom::error::ParseError<&'a [u8]>>> {
        let info = if self.original_unparsed_info == b"." || self.original_unparsed_info == b"" {
            (&b""[..], Vec::new())
        } else {
            separated_list1(tag(b";"), parse_info_tag)(self.original_unparsed_info)?
        };

        let input = self.other;
        let (input, format) = if !input.is_empty() {
            separated_list1(tag(b":"), take_while(|x| !is_separator2(x)))(input)?
        } else {
            (input, vec![])
        };

        let call_separator = |x: u8| match x {
            b',' | b':' | b'\t' | b'\n' => true,
            _ => false,
        };

        let (_input, call) = if !input.is_empty() {
            let (input, _) = tag(b"\t")(input)?; // consume tab
            separated_list1(
                tag(b"\t"),
                separated_list1(
                    tag(b":"),
                    separated_list1(tag(b","), take_till(call_separator)),
                ),
            )(input)?
        } else {
            (input, vec![])
        };

        Ok(CompleteVCFRecord {
            line: self.line,
            contig: self.contig,
            id: self.id,
            position: self.position,
            reference: self.reference,
            alternative: self.alternative,
            qual: self.qual,
            filter: self.filter,
            info: info
                .1
                .into_iter()
                .map(|InfoPair(x, y)| {
                    (Cow::Borrowed(x), y.into_iter().map(Cow::Borrowed).collect())
                })
                .collect(),
            format: format.into_iter().map(Cow::Borrowed).collect(),
            call: call
                .into_iter()
                .map(|x| {
                    x.into_iter()
                        .map(|y| y.into_iter().map(Cow::Borrowed).collect())
                        .collect()
                })
                .collect(),
        })
    }
}

fn is_separator(c: u8) -> bool {
    match c {
        b';' | b'\t' | b'\n' => true,
        _ => false,
    }
}

fn is_separator2(c: u8) -> bool {
    match c {
        b';' | b'\t' | b'\n' | b':' | b'=' | b',' => true,
        _ => false,
    }
}

fn parse_info_tag<'a>(input: &'a [u8]) -> IResult<&'a [u8], InfoPair> {
    let (input, info_data) = take_till(is_separator)(input)?;
    let (info_data, tag_name) = take_till(|x| x == b'=')(info_data)?;
    let data = if info_data.is_empty() {
        vec![]
    } else {
        let (info_data, _) = tag(b"=")(info_data)?;
        let (_info_data, data) = separated_list1(tag(b","), take_till(|x| x == b','))(info_data)?;
        data
    };

    Ok((input, InfoPair(tag_name, data)))
}

fn parse_vcf_internal<'a>(
    line: u32,
    input: &'a [u8],
) -> Result<PartialVCFRecord, nom::Err<(&'a [u8], nom::error::ErrorKind)>> {
    let raw = input;
    let (input, contig) = take_till(|c: u8| c == b'\t')(input)?;
    let (input, _) = tag(b"\t")(input)?; // consume tab
    let (input, position) = digit1(input)?;
    let (input, _) = tag(b"\t")(input)?; // consume tab
    let (input, id) = take_till(|c: u8| c == b'\t')(input)?;
    let (input, _) = tag(b"\t")(input)?; // consume tab
    let (input, reference) = take_till(|c: u8| c == b'\t')(input)?;
    let (input, _) = tag(b"\t")(input)?; // consume tab
    let (input, alternative) = separated_list1(
        tag(b","),
        take_while(|c: u8| match c {
            b',' | b'\t' => false,
            _ => true,
        }),
    )(input)?;

    // QUAL
    let (input, qual) = if !input.is_empty() && input != b"\n" {
        let (input, _) = tag(b"\t")(input)?; // consume tab
        take_till(|c: u8| c == b'\t' || c == b'\n')(input)?
    } else {
        (input, &[][..])
    };

    // FILTER
    let (input, filter) = if !input.is_empty() && input != b"\n" {
        let (input, _) = tag(b"\t")(input)?; // consume tab
        take_till(|c: u8| c == b'\t' || c == b'\n')(input)?
    } else {
        (input, &[][..])
    };

    // INFO
    let (input, unparsed_info) = if !input.is_empty() && input != b"\n" {
        let (input, _) = tag(b"\t")(input)?; // consume tab
        take_till(|c: u8| c == b'\t' || c == b'\n')(input)?
    } else {
        (input, &[][..])
    };

    // FORMAT/CALL
    let (input, other) = if !input.is_empty() && input != b"\n" {
        let (input, _) = tag(b"\t")(input)?; // consume tab
        take_while(|c: u8| c != b'\n')(input)?
    } else {
        (input, &[][..])
    };
    let (_input, _) = tag(b"\n")(input)?;

    Ok(PartialVCFRecord {
        raw,
        line,
        contig: Cow::Borrowed(contig),
        position: str::from_utf8(position).unwrap().parse().unwrap(),
        id: Cow::Borrowed(id),
        reference: Cow::Borrowed(reference),
        alternative: alternative.into_iter().map(Cow::Borrowed).collect(),
        qual: Cow::Borrowed(qual),
        filter: Cow::Borrowed(filter),
        original_unparsed_info: unparsed_info,
        unparsed_info: Cow::Borrowed(unparsed_info),
        other,
    })
}
