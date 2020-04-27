mod completerecord;
mod error;
mod partialrecord;
mod reader;
mod writer;

pub use completerecord::CompleteVCFRecord;
pub use error::{as_header_error, as_record_error, VCFParseError, VCFParseErrorKind};
pub use partialrecord::PartialVCFRecord;
pub use reader::VCFReader;
pub use writer::VCFWriter;

use nom::branch::alt;
use nom::bytes::complete::{tag, take_till, take_while};
use nom::multi::separated_list;
use nom::sequence::{delimited, separated_pair};
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt;
use std::io::{self, Write};
use std::str;

pub struct InfoPair<'a>(&'a [u8], Vec<&'a [u8]>);

pub trait VCFRecord: fmt::Debug + Clone + PartialEq {
    fn contig(&self) -> &[u8];
    fn position(&self) -> u64;
    fn id(&self) -> &[u8];
    fn reference(&self) -> &[u8];
    fn alternative(&self) -> &[Cow<[u8]>];
    fn qual(&self) -> &[u8];
    fn filter(&self) -> &[u8];
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()>;
}

#[derive(PartialEq, Clone)]
pub struct VCFHeaderItem {
    pub key: Vec<u8>,
    pub value: Vec<u8>,
    pub detail: HashMap<Vec<u8>, Vec<u8>>,
}

impl fmt::Debug for VCFHeaderItem {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(f, "VCFHeaderItem {{ ")?;
        write!(f, "key: b\"{}\", ", str::from_utf8(&self.key).unwrap())?;
        write!(f, "value: b\"{}\", ", str::from_utf8(&self.value).unwrap())?;
        write!(f, "detail: {{ ")?;
        for (i, (k, v)) in self.detail.iter().enumerate() {
            if i != 0 {
                write!(f, ", ")?;
            }
            write!(
                f,
                "b\"{}\": b\"{}\"",
                str::from_utf8(k).unwrap(),
                str::from_utf8(v).unwrap()
            )?;
        }
        write!(f, " }}")?;
        write!(f, " }}")?;
        Ok(())
    }
}

impl VCFHeaderItem {
    pub fn parse(input: &[u8], line: u32) -> Result<VCFHeaderItem, VCFParseError> {
        match VCFHeaderItem::parse_helper(input) {
            Ok(item) => Ok(item),
            Err(nom::Err::Failure(e)) | Err(nom::Err::Error(e)) => {
                Err(as_header_error(line, input, e))
            }
            _ => unreachable!(),
        }
    }

    fn parse_helper(
        input: &[u8],
    ) -> Result<VCFHeaderItem, nom::Err<(&[u8], nom::error::ErrorKind)>> {
        let (input, _) = tag(b"##")(input)?;
        let (input, key) = take_till(|x| x == b'=')(input)?;
        let (input, _) = tag(b"=")(input)?;
        let (_input, value) = take_till(|x| x == b'\n')(input)?;
        let mut detail = HashMap::new();
        if value.starts_with(b"<") {
            let (_, detail_vec) = delimited(
                tag("<"),
                separated_list(
                    tag(b","),
                    separated_pair(
                        take_while(|x| match x {
                            b'=' | b',' | b'>' => false,
                            _ => true,
                        }),
                        tag(b"="),
                        alt((
                            delimited(
                                tag("\""),
                                take_while(|x| match x {
                                    b'"' => false,
                                    _ => true,
                                }),
                                tag("\""),
                            ),
                            take_while(|x| match x {
                                b'=' | b',' | b'>' => false,
                                _ => true,
                            }),
                        )),
                    ),
                ),
                tag(">"),
            )(value)?;
            for one in detail_vec {
                detail.insert(one.0.to_vec(), one.1.to_vec());
            }
        }

        Ok(VCFHeaderItem {
            key: key.to_vec(),
            value: value.to_vec(),
            detail,
        })
    }

    pub fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(b"##")?;
        writer.write_all(&self.key)?;
        writer.write_all(b"=")?;
        writer.write_all(&self.value)?;
        writer.write_all(b"\n")?;
        Ok(())
    }
}

#[derive(PartialEq, Clone)]
pub struct VCFHeader {
    pub header_items: Vec<VCFHeaderItem>,
    pub samples: Vec<Vec<u8>>,
}

#[cfg(test)]
mod test;
