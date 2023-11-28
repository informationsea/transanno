use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<impl Read> {
    Ok(autocompress::autodetect_open(path)?)
}

pub fn create<P: AsRef<Path>>(path: P) -> std::io::Result<impl Write> {
    let extension = path.as_ref().extension().and_then(|x| x.to_str());
    let writer: Box<dyn Write> = match extension {
        Some("gz") => Box::new(bgzip::BGZFWriter::new(
            File::create(path)?,
            bgzip::Compression::default(),
        )),
        _ => Box::new(autocompress::autodetect_create(
            path,
            autocompress::CompressionLevel::default(),
        )?),
    };
    Ok(writer)
}
