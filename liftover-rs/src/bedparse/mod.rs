use nom::{
    branch::alt,
    bytes::complete::{tag, take_till1, take_while1},
    character::is_digit,
    combinator::{eof, map, map_res, opt},
    sequence::tuple,
    IResult,
};
use std::io::{self, prelude::*};
use std::str;

#[derive(Debug, PartialEq, Eq)]
pub struct BedEntry<'a> {
    pub chrom: &'a [u8],
    pub start: u64,
    pub end: u64,
    pub name: Option<&'a [u8]>,
    pub score: Option<&'a [u8]>,
    pub strand: Option<&'a [u8]>,
    pub others: Option<&'a [u8]>,
}

impl<'a> BedEntry<'a> {
    pub fn new(
        chrom: &'a [u8],
        start: u64,
        end: u64,
        name: Option<&'a [u8]>,
        score: Option<&'a [u8]>,
        strand: Option<&'a [u8]>,
        others: Option<&'a [u8]>,
    ) -> Self {
        Self {
            chrom,
            start,
            end,
            name,
            score,
            strand,
            others,
        }
    }

    pub fn with_new_name(&self, new_name: &'a [u8]) -> Self {
        Self {
            chrom: self.chrom,
            start: self.start,
            end: self.end,
            name: Some(new_name),
            score: self.score,
            strand: self.strand,
            others: self.others,
        }
    }

    pub fn with_new_position(
        &self,
        chrom: &'a [u8],
        start: u64,
        end: u64,
        strand: Option<&'a [u8]>,
    ) -> Self {
        Self {
            chrom,
            start,
            end,
            name: self.name,
            score: self.score,
            strand,
            others: self.others,
        }
    }

    pub fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(self.chrom)?;
        writer.write_all(b"\t")?;
        writer.write_all(self.start.to_string().as_bytes())?;
        writer.write_all(b"\t")?;
        writer.write_all(self.end.to_string().as_bytes())?;
        if let Some(name) = self.name {
            writer.write_all(b"\t")?;
            writer.write_all(name)?;
            if let Some(score) = self.score {
                writer.write_all(b"\t")?;
                writer.write_all(score)?;
                if let Some(strand) = self.strand {
                    writer.write_all(b"\t")?;
                    writer.write_all(strand)?;
                    if let Some(others) = self.others {
                        writer.write_all(b"\t")?;
                        writer.write_all(others)?;
                    }
                }
            }
        }
        writer.write_all(b"\n")?;
        Ok(())
    }
}

pub fn opt_field<'a, O, E, F>(
    opt_parser: F,
) -> impl FnMut(&'a [u8]) -> IResult<&'a [u8], Option<(&'a [u8], O)>, E>
where
    F: nom::Parser<&'a [u8], O, E>,
    E: nom::error::ParseError<&'a [u8]>,
{
    opt(map(
        tuple((
            tag(b"\t"),
            take_while1(|x| x != b'\t' && x != b'\n'),
            opt_parser,
        )),
        |(_, x, y)| (x, y),
    ))
}

pub fn parse_bed_line<'a>(input: &'a [u8]) -> IResult<&'a [u8], BedEntry<'a>> {
    let (input, chrom) = take_till1(|x| x == b'\t')(input)?;
    let (input, _) = tag(b"\t")(input)?;
    let (input, start) = map_res(
        map_res(take_while1(is_digit), str::from_utf8),
        str::parse::<u64>,
    )(input)?;
    let (input, _) = tag(b"\t")(input)?;
    let (input, end) = map_res(
        map_res(take_while1(is_digit), str::from_utf8),
        str::parse::<u64>,
    )(input)?;

    let (input, (others, _)) = tuple((
        opt_field(opt_field(opt_field(opt(tuple((
            tag(b"\t"),
            take_while1(|x| x != b'\n'),
        )))))),
        alt((eof, tag("\n"))),
    ))(input)?;

    Ok((
        input,
        BedEntry {
            chrom,
            start,
            end,
            name: others.map(|x| x.0),
            score: others.map(|x| x.1.map(|x| x.0)).flatten(),
            strand: others
                .map(|x| x.1.map(|x| x.1.map(|x| x.0)))
                .flatten()
                .flatten(),
            others: others
                .map(|x| x.1.map(|x| x.1.map(|x| x.1.map(|x| x.1))))
                .flatten()
                .flatten()
                .flatten(),
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse1() {
        let input = b"chr1\t1\t2\t3\n";
        let (_, entry) = parse_bed_line(input).unwrap();
        assert_eq!(
            entry,
            BedEntry {
                chrom: b"chr1",
                start: 1,
                end: 2,
                name: Some(b"3"),
                score: None,
                strand: None,
                others: None,
            }
        );
        let mut output = Vec::new();
        entry.write(&mut output).unwrap();
        assert_eq!(output, input);
    }

    #[test]
    fn test_parse2() {
        let input = b"chr1\t1\t2";
        let (_, entry) = parse_bed_line(input).unwrap();
        assert_eq!(
            entry,
            BedEntry {
                chrom: b"chr1",
                start: 1,
                end: 2,
                name: None,
                score: None,
                strand: None,
                others: None,
            }
        );
        let mut output = Vec::new();
        entry.write(&mut output).unwrap();
        assert_eq!(output, b"chr1\t1\t2\n");
    }

    #[test]
    fn test_parse3() {
        let input = b"chr1\t1\t2\n";
        let (_, entry) = parse_bed_line(input).unwrap();
        assert_eq!(
            entry,
            BedEntry {
                chrom: b"chr1",
                start: 1,
                end: 2,
                name: None,
                score: None,
                strand: None,
                others: None,
            }
        );
        let mut output = Vec::new();
        entry.write(&mut output).unwrap();
        assert_eq!(output, input);
    }

    #[test]
    fn test_parse4() {
        let input = b"chr1\t1\t2\nx";
        let (remain, entry) = parse_bed_line(input).unwrap();
        assert_eq!(
            entry,
            BedEntry {
                chrom: b"chr1",
                start: 1,
                end: 2,
                name: None,
                score: None,
                strand: None,
                others: None,
            }
        );
        assert_eq!(remain, b"x");
        let mut output = Vec::new();
        entry.write(&mut output).unwrap();
        assert_eq!(output, b"chr1\t1\t2\n");
    }

    #[test]
    fn test_parse6() {
        let input = b"chr1\t1\t2\tname\tscore\tstrand\tothers1\tothers2\n";
        let (remain, entry) = parse_bed_line(input).unwrap();
        assert_eq!(
            entry,
            BedEntry {
                chrom: b"chr1",
                start: 1,
                end: 2,
                name: Some(b"name"),
                score: Some(b"score"),
                strand: Some(b"strand"),
                others: Some(b"others1\tothers2"),
            }
        );
        assert_eq!(remain, b"");
        let mut output = Vec::new();
        entry.write(&mut output).unwrap();
        assert_eq!(output, input);
    }

    #[test]
    fn test_parse7() {
        let input = b"chr1\t1\t2\tname\tscore\tstrand\n";
        let (remain, entry) = parse_bed_line(input).unwrap();
        assert_eq!(
            entry,
            BedEntry {
                chrom: b"chr1",
                start: 1,
                end: 2,
                name: Some(b"name"),
                score: Some(b"score"),
                strand: Some(b"strand"),
                others: None,
            }
        );
        assert_eq!(remain, b"");
        let mut output = Vec::new();
        entry.write(&mut output).unwrap();
        assert_eq!(output, input);
    }

    #[test]
    fn test_parse8() {
        let input = b"chr1\t1\t2\tname\tscore\n";
        let (remain, entry) = parse_bed_line(input).unwrap();
        assert_eq!(
            entry,
            BedEntry {
                chrom: b"chr1",
                start: 1,
                end: 2,
                name: Some(b"name"),
                score: Some(b"score"),
                strand: None,
                others: None,
            }
        );
        assert_eq!(remain, b"");
        let mut output = Vec::new();
        entry.write(&mut output).unwrap();
        assert_eq!(output, input);
    }
}
