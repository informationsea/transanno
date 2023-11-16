use std::io::{Read, Write};
use std::path::Path;

use autocompress::io::{RayonReader, RayonWriter};

pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<impl Read> {
    Ok(RayonReader::new(autocompress::autodetect_open(path)?))
}

pub fn create<P: AsRef<Path>>(path: P) -> std::io::Result<impl Write> {
    Ok(RayonWriter::new(
        autocompress::autodetect_create_prefer_bgzip(
            path,
            autocompress::CompressionLevel::default(),
        )?,
    ))
}
