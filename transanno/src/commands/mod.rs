mod chain2bedvcf;
// mod chain2chunkbed;
mod chain_left_align;
// mod generate_completions;
mod liftbed;
mod liftgene;
mod liftvcf;
mod minimap2chain;

#[derive(Debug, Clone, clap::Subcommand)]
pub enum Commands {
    LeftAlign(chain_left_align::ChainLeftAlign),
    ChainToBedVcf(chain2bedvcf::Chain2BedVcf),
    Liftgene(liftgene::LiftGene),
    Minimap2chain(minimap2chain::Minimap2Chain),
    Liftvcf(liftvcf::LiftVcf),
    Liftbed(liftbed::LiftBed),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::LeftAlign(x) => x.run(),
            Commands::ChainToBedVcf(x) => x.run(),
            Commands::Liftgene(x) => x.run(),
            Commands::Minimap2chain(x) => x.run(),
            Commands::Liftvcf(x) => x.run(),
            Commands::Liftbed(x) => x.run(),
        }
    }
}
