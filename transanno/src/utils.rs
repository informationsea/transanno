use std::io::{BufRead, BufReader, Write};
use std::path::Path;

pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<Box<dyn BufRead>> {
    let path = path.as_ref();
    let mut file = BufReader::new(std::fs::File::open(path)?);
    let buf = file.fill_buf()?;
    if buf[0] == 0x00 && buf[1] == 0x00 {
        Ok(Box::new(BufReader::new(
            flate2::bufread::MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(file))
    }
}

pub fn create<P: AsRef<Path>>(path: P) -> std::io::Result<Box<dyn Write>> {
    let path = path.as_ref();
    let file = std::fs::File::create(path)?;
    if let Some(extension) = path.extension() {
        if extension == "gz" {
            return Ok(Box::new(bgzip::BGZFWriter::new(
                file,
                bgzip::Compression::default(),
            )));
        }
    }
    Ok(Box::new(file))
}
