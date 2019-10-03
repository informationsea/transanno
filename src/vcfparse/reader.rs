use super::{
    as_header_error, PartialVCFRecord, VCFHeader, VCFHeaderItem, VCFParseError, VCFParseErrorKind,
};
use nom::bytes::complete::{tag, take_till};
use nom::multi::many0;
use std::io::{self, BufRead};

pub struct VCFReader<R: BufRead> {
    reader: R,
    pub header: VCFHeader,
    current_line: u32,
    buffer: Vec<u8>,
}

impl<R: io::Read> VCFReader<io::BufReader<R>> {
    pub fn new(reader: R) -> Result<VCFReader<io::BufReader<R>>, VCFParseError> {
        let mut reader = io::BufReader::new(reader);
        let mut current_line = 0;
        let mut header_items = Vec::new();
        let mut samples = Vec::new();
        let mut line = Vec::new();
        loop {
            line.clear();
            reader.read_until(b'\n', &mut line)?;
            current_line += 1;

            if line.starts_with(b"##") {
                header_items.push(VCFHeaderItem::parse(&line, current_line)?);
            } else if line.starts_with(b"#") {
                let parsed_samples = parse_samples(&line).map_err(|e| match e {
                    nom::Err::Error(e) | nom::Err::Failure(e) => {
                        as_header_error(current_line, &line, e)
                    }
                    _ => unreachable!(),
                })?;
                for one in parsed_samples {
                    samples.push(one.to_vec());
                }
                break;
            } else {
                return Err(VCFParseErrorKind::HeaderParseError {
                    line: current_line,
                    column: 0,
                }
                .into());
            }
        }

        Ok(VCFReader {
            reader,
            header: VCFHeader {
                header_items,
                samples,
            },
            current_line,
            buffer: Vec::new(),
        })
    }

    pub fn next_record(&mut self) -> Result<Option<PartialVCFRecord>, VCFParseError> {
        self.buffer.clear();
        self.current_line += 1;
        let read_bytes = self.reader.read_until(b'\n', &mut self.buffer)?;
        if read_bytes == 0 {
            Ok(None)
        } else {
            PartialVCFRecord::parse_vcf(self.current_line, &self.buffer).map(Some)
        }
    }
}

type NomU8Err<'a> = nom::Err<(&'a [u8], nom::error::ErrorKind)>;
fn parse_samples(input: &[u8]) -> Result<Vec<&[u8]>, NomU8Err> {
    let (input, _) = tag("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")(input)?;
    let samples = Vec::new();
    if input.is_empty() || input == b"\n" || input == b"\r" {
        return Ok(samples);
    }
    let (input, _) = tag("\t")(input)?;
    if input.is_empty() || input == b"\n" || input == b"\r" {
        return Ok(samples);
    }
    let (input, _) = tag("FORMAT")(input)?;
    if input.is_empty() || input == b"\n" || input == b"\r" {
        return Ok(samples);
    }
    let (_input, samples) = many0(|x| {
        let (x, _) = tag(b"\t")(x)?;
        take_till(|y| match y {
            b'\t' | b'\n' | b'\r' => true,
            _ => false,
        })(x)
    })(input)?;

    Ok(samples)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_samples() {
        assert_eq!(
            parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n").unwrap(),
            Vec::<&[u8]>::new()
        );
        assert_eq!(
            parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n").unwrap(),
            Vec::<&[u8]>::new()
        );

        assert_eq!(
            parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfoo\n").unwrap(),
            vec![b"foo"]
        );

        assert_eq!(
            parse_samples(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfoo\tbar\n")
                .unwrap(),
            vec![b"foo", b"bar"]
        );
    }
}
