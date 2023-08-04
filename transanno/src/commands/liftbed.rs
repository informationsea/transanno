use anyhow::Context;
use autocompress::{create, open, CompressionLevel};
use clap::Args;
use liftover::{bedparse::parse_bed_line, chain::Strand, poslift::PositionLiftOver};
use std::io::{self, prelude::*};

#[derive(Debug, Clone, Args)]
#[command(about = "Lift BED file")]
pub struct LiftBed {
    #[arg(long = "chain", short = 'c', help = "chain file")]
    chain: String,
    #[arg(help = "input BED file")]
    bed: String,
    #[arg(long = "output", short = 'o', help = "BED output path (unsorted)")]
    output: String,
    #[arg(long = "failed", short = 'f', help = "Failed BED output path")]
    failed: Option<String>,
    #[arg(help = "Allow multi-map", long = "allow-multi-map", short = 'm')]
    allow_multimap: bool,
}

impl LiftBed {
    pub fn run(&self) -> anyhow::Result<()> {
        let chain_file =
            PositionLiftOver::load(open(&self.chain).context("Failed to open chain file")?)
                .context("Failed parse chain file.")?;
        let mut writer = io::BufWriter::new(
            create(&self.output, CompressionLevel::Default)
                .with_context(|| format!("Failed to create {}", &self.output))?,
        );
        let mut failed_writer = if let Some(x) = &self.failed {
            Some(io::BufWriter::new(
                create(x, CompressionLevel::Default)
                    .with_context(|| format!("Failed to create {}", x))?,
            ))
        } else {
            None
        };

        let mut reader = io::BufReader::new(
            open(&self.bed).with_context(|| format!("Failed to open {}", &self.bed))?,
        );

        let mut line = Vec::new();
        let mut count = 0;
        while reader.read_until(b'\n', &mut line)? > 0 {
            count += 1;
            if line.starts_with(b"#") {
                writer.write_all(&line)?;
                line.clear();
                continue;
            }

            let (remain, bed) = parse_bed_line(&line)
                .map_err(|e| anyhow::anyhow!("Failed to parse BED at line: {} / {}", count, e))?;
            if !remain.is_empty() {
                return Err(anyhow::anyhow!(
                    "Failed to parse BED at line: {} / {}",
                    count,
                    String::from_utf8_lossy(&line)
                ));
            }

            let lift_results = chain_file.lift_region(
                String::from_utf8_lossy(bed.chrom).as_ref(),
                bed.start..bed.end,
            );
            if lift_results.is_empty() {
                if let Some(failed_writer) = failed_writer.as_mut() {
                    let new_name = format!(
                        "{};no-map",
                        bed.name
                            .map(|x| String::from_utf8_lossy(x).to_string())
                            .unwrap_or("".to_string())
                    );
                    let new_bed = bed.with_new_name(new_name.as_bytes());
                    new_bed.write(failed_writer)?;
                }
            } else if lift_results.len() > 1 && !self.allow_multimap {
                if let Some(failed_writer) = failed_writer.as_mut() {
                    let new_name = format!(
                        "{};multi-map",
                        bed.name
                            .map(|x| String::from_utf8_lossy(x).to_string())
                            .unwrap_or("".to_string())
                    );
                    let new_bed = bed.with_new_name(new_name.as_bytes());
                    new_bed.write(failed_writer)?;
                }
            } else {
                for one in lift_results {
                    let strand = match bed.strand {
                        Some(b"+") => match one.strand {
                            Strand::Forward => Some(&b"+"[..]),
                            Strand::Reverse => Some(&b"-"[..]),
                        },
                        Some(b"-") => match one.strand {
                            Strand::Forward => Some(&b"-"[..]),
                            Strand::Reverse => Some(&b"+"[..]),
                        },
                        _ => bed.strand,
                    };
                    let new_bed = bed.with_new_position(
                        one.chromosome.name.as_bytes(),
                        one.start,
                        one.end,
                        strand,
                    );
                    new_bed.write(&mut writer)?;
                }
            }

            line.clear();
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::Cli;
    use clap::Parser;

    #[test]
    fn test_liftbed() {
        let cli = Cli::parse_from(&[
            "transanno",
            "liftbed",
            "--chain",
            "testfiles/GRCh38.primary_assembly-to-GRCh37.primary_assembly.chain.gz",
            "--output",
            "../target/mapped-to-hg19.bed",
            "--failed",
            "../target/failed-to-map-hg19.bed",
            "testfiles/hg38.bed",
        ]);
        cli.command.run().unwrap();
    }
}
