use super::VCFRecord;
use std::borrow::Cow;
use std::fmt::{self, Debug};
use std::io;
use std::io::Write;
use std::str;

type InfoCowPair<'a> = (Cow<'a, [u8]>, Vec<Cow<'a, [u8]>>);

#[derive(PartialEq, Clone)]
pub struct CompleteVCFRecord<'a> {
    pub line: u32,
    pub contig: Cow<'a, [u8]>,
    pub position: u64,
    pub id: Cow<'a, [u8]>,
    pub reference: Cow<'a, [u8]>,
    pub alternative: Vec<Cow<'a, [u8]>>,
    pub qual: Cow<'a, [u8]>,
    pub filter: Cow<'a, [u8]>,
    pub info: Vec<InfoCowPair<'a>>,
    pub format: Vec<Cow<'a, [u8]>>,
    pub call: Vec<Vec<Vec<Cow<'a, [u8]>>>>,
}

impl<'a> VCFRecord for CompleteVCFRecord<'a> {
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
        writer.write_all(&self.qual)?;
        writer.write_all(b"\t")?;
        writer.write_all(&self.filter)?;

        let should_write_call_format = !self.call.is_empty() || !self.format.is_empty();
        let should_write_info = !self.info.is_empty() || should_write_call_format;

        if should_write_info {
            writer.write_all(b"\t")?;
            for (i, one) in self.info.iter().enumerate() {
                if i != 0 {
                    writer.write_all(b";")?;
                }
                writer.write_all(&one.0)?;
                if !one.1.is_empty() {
                    writer.write_all(b"=")?;
                    for (j, y) in one.1.iter().enumerate() {
                        if j != 0 {
                            writer.write_all(b",")?;
                        }
                        writer.write_all(y)?;
                    }
                }
            }
        }

        if should_write_call_format {
            writer.write_all(b"\t")?;
            for (i, one) in self.format.iter().enumerate() {
                if i != 0 {
                    writer.write_all(b":")?;
                }
                writer.write_all(one)?;
            }
            for one_sample in self.call.iter() {
                writer.write_all(b"\t")?;
                for (i, one_call) in one_sample.iter().enumerate() {
                    if i != 0 {
                        writer.write_all(b":")?;
                    }
                    for (j, one_value) in one_call.iter().enumerate() {
                        if j != 0 {
                            writer.write_all(b",")?;
                        }
                        writer.write_all(one_value)?;
                    }
                }
            }
        }
        writer.write_all(b"\n")?;
        Ok(())
    }
}

impl<'a> Debug for CompleteVCFRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        write!(f, "CompleteVCFRecord {{ ")?;
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

        write!(f, "info: vec![")?;
        for (i, x) in self.info.iter().enumerate() {
            if i != 0 {
                write!(f, ", ")?;
            }
            write!(f, "(b\"{}\", vec![", str::from_utf8(&x.0).unwrap())?;
            for (j, y) in x.1.iter().enumerate() {
                if j != 0 {
                    write!(f, ", ")?;
                }
                write!(f, "b\"{}\"", str::from_utf8(y).unwrap())?;
            }
            write!(f, "])")?;
        }
        write!(f, "], ")?;

        write!(f, "format: vec![")?;
        for (i, x) in self.format.iter().enumerate() {
            if i != 0 {
                write!(f, ", ")?;
            }
            write!(f, "b\"{}\"", str::from_utf8(x).unwrap())?;
        }
        write!(f, "], ")?;

        write!(f, "call: vec![")?;
        for (i, x) in self.call.iter().enumerate() {
            if i != 0 {
                write!(f, ", ")?;
            }
            write!(f, "vec![")?;
            for (j, y) in x.iter().enumerate() {
                if j != 0 {
                    write!(f, ", ")?;
                }
                write!(f, "vec![")?;
                for (k, z) in y.iter().enumerate() {
                    if k != 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "b\"{}\"", str::from_utf8(z).unwrap())?;
                }
                write!(f, "]")?;
            }
            write!(f, "]")?;
        }
        write!(f, "] }}")?;
        Ok(())
    }
}
