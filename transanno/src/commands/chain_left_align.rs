use anyhow::Context;
use autocompress::{create, open, CompressionLevel};
use bio::io::fasta::IndexedReader;
use clap::Args;
use liftover::chain;
use log::info;

#[derive(Debug, Clone, Args)]
#[command(about = "Left align and normalize chain file")]
pub struct ChainLeftAlign {
    #[arg(help = "Original chain file")]
    original_chain: String,
    #[arg(help = "Output chain file", short = 'o', long = "output")]
    output: String,
    #[arg(
        help = "Original assembly FASTA (.fai file is required)",
        short = 'r',
        long = "original",
        alias = "reference"
    )]
    original_sequence: String,
    #[arg(
        help = "New assembly FASTA (.fai file is required)",
        short = 'q',
        long = "new",
        alias = "query"
    )]
    new_sequence: String,
}

impl ChainLeftAlign {
    pub fn run(&self) -> anyhow::Result<()> {
        info!("start loading chain");
        let chain_file = open(&self.original_chain).context("Cannot create input chain file")?;
        let mut output_file = create(&self.output, CompressionLevel::Default)
            .context("Cannot create output chain file")?;
        let mut original_seq = IndexedReader::from_file(&self.original_sequence)
            .context("Cannot load original assembly FASTA")?;
        let mut new_seq = IndexedReader::from_file(&self.new_sequence)
            .context("Cannot load new assembly FASTA")?;
        let chain_data = chain::ChainFile::load(chain_file).expect("Failed to parse chain file");
        let left_aligned = chain_data.left_align(&mut original_seq, &mut new_seq)?;
        left_aligned.write(&mut output_file)?;
        Ok(())
    }
}
