mod chain2bedvcf;
// mod chain2chunkbed;
mod chain_left_align;
// mod generate_completions;
mod liftgene;
mod liftvcf;
mod minimap2chain;

// pub use chain2bedvcf::chain_to_bed_vcf;
// pub use liftgene::lift_gene;
// pub use minimap2chain::minimap2_to_chain;

// use clap::{crate_authors, crate_version, App, ArgMatches, SubCommand};

// pub(crate) const COMMANDS: &[&dyn Command] = &[
//     &minimap2chain::Minimap2Chain,
//     &liftgene::LiftGene,
//     &liftvcf::LiftVcf,
//     &chain2bedvcf::Chain2BedVcf,
//     //&chain2chunkbed::Chain2ChunkBed, // Developing
//     &chain_left_align::ChainLeftAlign,
//     &generate_completions::GenerateCompletions,
// ];

// pub trait Command {
//     fn cli(&self) -> App<'static, 'static> {
//         self.config_subcommand(SubCommand::with_name(self.command_name()))
//             .version(crate_version!())
//             .author(crate_authors!())
//     }
//     fn command_name(&self) -> &'static str;
//     fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static>;
//     fn run(&self, matches: &ArgMatches<'static>) -> anyhow::Result<()>;
// }

#[derive(Debug, Clone, clap::Subcommand)]
pub enum Commands {
    LeftAlign(chain_left_align::ChainLeftAlign),
    ChainToBedVcf(chain2bedvcf::Chain2BedVcf),
    Liftgene(liftgene::LiftGene),
    Minimap2chain(minimap2chain::Minimap2Chain),
    Liftvcf(liftvcf::LiftVcf),
}

impl Commands {
    pub fn run(&self) -> anyhow::Result<()> {
        match self {
            Commands::LeftAlign(x) => x.run(),
            Commands::ChainToBedVcf(x) => x.run(),
            Commands::Liftgene(x) => x.run(),
            Commands::Minimap2chain(x) => x.run(),
            Commands::Liftvcf(x) => x.run(),
        }
    }
}
