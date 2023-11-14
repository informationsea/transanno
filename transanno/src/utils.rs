use std::io::{Read, Write};
use std::path::Path;

pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<impl Read> {
    autocompress::autodetect_open(path)
}

pub fn create<P: AsRef<Path>>(path: P) -> std::io::Result<impl Write> {
    autocompress::autodetect_create_prefer_bgzip(path, autocompress::CompressionLevel::default())
}
