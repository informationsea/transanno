use super::{VCFHeader, VCFParseError, VCFRecord};
use std::io;
use std::io::prelude::*;

pub struct VCFWriter<W: Write> {
    pub header: VCFHeader,
    writer: W,
    buffer: Vec<u8>,
}

impl<W: Write> VCFWriter<io::BufWriter<W>> {
    pub fn new(writer: W, header: VCFHeader) -> Result<VCFWriter<io::BufWriter<W>>, VCFParseError> {
        let mut writer = io::BufWriter::new(writer);
        for one_item in header.header_items.iter() {
            one_item.write(&mut writer)?;
        }
        if header.samples.is_empty() {
            writer.write_all(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
        } else {
            writer.write_all(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
            for one in header.samples.iter() {
                writer.write_all(b"\t")?;
                writer.write_all(one)?;
            }
        }
        writer.write_all(b"\n")?;

        Ok(VCFWriter {
            header,
            writer,
            buffer: Vec::new(),
        })
    }
}

impl<W: Write> VCFWriter<W> {
    pub fn write_record<R: VCFRecord>(&mut self, record: &R) -> Result<(), VCFParseError> {
        self.buffer.clear();
        record.write(&mut self.buffer)?;
        self.writer.write_all(&self.buffer)?;
        Ok(())
    }
}
