use anyhow::Context;
use autocompress::{create, open, CompressionLevel};
use bio::io::fasta::IndexedReader;
use clap::Args;
use liftover::{chain, variantlift, vcflift, LiftOverError};
use log::info;

#[derive(Debug, Clone, Args)]
#[command(about = "LiftOver VCF file")]
pub struct LiftVcf {
    #[arg(
        long = "original-assembly",
        alias = "reference",
        short = 'r',
        help = "Original assembly FASTA (.fai file is required)"
    )]
    original_sequence: String,
    #[arg(
        long = "new-assembly",
        alias = "query",
        short = 'q',
        help = "New assembly FASTA (.fai file is required)"
    )]
    new_sequence: String,
    #[arg(long = "chain", short = 'c', help = "chain file")]
    chain: String,
    #[arg(long = "vcf", short = 'v', help = "input VCF file to liftOver")]
    vcf: String,
    #[arg(
        long = "output",
        short = 'o',
        help = "output VCF file for succeeded to liftOver records (This file is not sorted)"
    )]
    output: String,
    #[arg(
        long = "fail",
        short = 'f',
        help = "output VCF file for failed to liftOver records"
    )]
    fail: String,
    #[arg(help = "Allow multi-map", long = "allow-multi-map", short = 'm')]
    allow_multimap: bool,
    #[arg(
        help = "length of acceptable deletion",
        long = "acceptable-deletion",
        short = 'd',
        default_value = "3"
    )]
    acceptable_deletion: u64,
    #[arg(
        help = "length of acceptable insertion",
        long = "acceptable-insertion",
        short = 'i',
        default_value = "3"
    )]
    acceptable_insertion: u64,
    #[arg(
        long = "no-rewrite-format",
        help = "Do not rewrite order of FORMAT tags"
    )]
    do_not_rewrite_info: bool,
    #[arg(long = "no-rewrite-gt", help = "Do not rewrite order of GT")]
    do_not_rewrite_gt: bool,
    #[arg(
        long = "no-rewrite-allele-frequency",
        help = "Do not rewrite AF or other allele frequency info"
    )]
    do_not_rewrite_allele_frequency: bool,
    #[arg(
        long = "no-rewrite-allele-count",
        help = "Do not rewrite AC or other count frequency info"
    )]
    do_not_rewrite_allele_count: bool,
    #[arg(
        long = "noswap",
        help = "Do not swap ref/alt when reference allele is changed. This option is suitable to do liftOver clinVar, COSMIC annotations"
    )]
    do_not_swap_ref_alt: bool,
    #[arg(
        long = "no-left-align-chain",
        help = "Do not run left align chain file"
    )]
    do_not_left_align_chain: bool,
    #[arg(
        long = "do-not-use-dot-when-alt-equal-to-ref",
        help = "Do not use dot as ALT when ALT column is equal to REF"
    )]
    do_not_use_dot_when_alt_equal_to_ref: bool,
    #[arg(
        long = "do-not-prefer-same-contig-when-multimap",
        help = "Do not prefer same name contig when a variant lifted into multiple positions. (When you use this option, a variant which lifted into a main chromosome and alternative contigs, lift over will be failed if multimap is not allowed)"
    )]
    do_not_prefer_cis_contig_when_multimap: bool,
    #[arg(
        long = "ignore-fasta-length-mismatch",
        help = "Ignore length mismatch between chain and fasta file"
    )]
    ignore_fasta_length_mismatch: bool,
}

impl LiftVcf {
    pub fn run(&self) -> anyhow::Result<()> {
        info!("start loading chain and fasta");
        let mut original_seq = IndexedReader::from_file(&self.original_sequence)
            .context("Failed to load original assembly FASTA")?;
        let mut new_seq = IndexedReader::from_file(&self.new_sequence)
            .context("Failed to load new assembly FASTA")?;
        let chain = chain::ChainFile::load(autocompress::open(&self.chain)?)?
            .left_align(&mut original_seq, &mut new_seq)
            .context("Failed to load chain file")?;
        // Reference/Query sequence and chain consistency
        for one_chain in chain.chain_list.iter() {
            match one_chain.check_sequence_consistency(&mut original_seq, &mut new_seq) {
                Ok(_) => (),
                Err(e) => match e {
                    LiftOverError::ChromosomeNotFound(_) => (),
                    LiftOverError::UnmatchedChromosomeLength(_, _, _)
                    | LiftOverError::QueryChromosomeLengthIsNotMatch(_)
                    | LiftOverError::ReferenceChromosomeLengthIsNotMatch(_) => {
                        if !self.ignore_fasta_length_mismatch {
                            return Err(e.into());
                        }
                    }
                    _ => return Err(e.into()),
                },
            }
        }

        let variant_liftover = variantlift::VariantLiftOver::new(chain, original_seq, new_seq);
        let mut vcf_lift = vcflift::VCFLiftOver::new(
            variant_liftover,
            vcflift::VCFLiftOverParameters::new()
                .allow_multimap(self.allow_multimap)
                .acceptable_deletion(self.acceptable_deletion)
                .acceptable_insertion(self.acceptable_insertion)
                .do_not_rewrite_info(self.do_not_rewrite_info)
                //.do_not_rewrite_format(self.do_not_)
                .do_not_rewrite_gt(self.do_not_rewrite_gt)
                .do_not_rewrite_allele_frequency(self.do_not_rewrite_allele_frequency)
                .do_not_rewrite_allele_count(self.do_not_rewrite_allele_count)
                .do_not_swap_ref_alt(self.do_not_swap_ref_alt)
                .do_not_left_align_chain_file(self.do_not_left_align_chain)
                .do_not_use_dot_when_alt_equal_to_ref(self.do_not_use_dot_when_alt_equal_to_ref)
                .do_not_prefer_cis_contig_when_multimap(
                    self.do_not_prefer_cis_contig_when_multimap,
                ),
        );
        info!("chain file and fasta files were loaded");

        let uncompressed_reader = open(&self.vcf).expect("Cannot open input VCF");
        let success_writer =
            create(&self.output, CompressionLevel::Default).expect("Cannot create output VCF");
        let failed_writer = create(&self.fail, CompressionLevel::Default)
            .expect("Cannot create output VCF for failed records");
        vcf_lift.lift_vcf(uncompressed_reader, success_writer, failed_writer)?;
        Ok(())
    }
}
