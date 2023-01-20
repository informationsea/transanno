use anyhow::Context;
use autocompress::{create, open, CompressionLevel};
use clap::Args;
use once_cell::sync::Lazy;
use regex::Regex;
use std::io::{self, BufRead, Write};

#[derive(Debug, Clone, Args)]
#[command(
    about = "Convert minimap2 result to chain file",
    long_about = r#"Convert minimap2 result to chain file

A paf file should be created with a command shown in below.

$ minimap2 -cx asm5 --cs NEW_FASTA ORIGINAL_FASTA > PAF_FILE.paf
"#
)]
pub struct Minimap2Chain {
    #[arg(help = "Input paf file")]
    paf: String,
    #[arg(long = "output", short = 'o', help = "Output chain file")]
    output: String,
}

impl Minimap2Chain {
    //     fn command_name(&self) -> &'static str {
    //         "minimap2chain"
    //     }
    //     fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
    //         app.about("Convert minimap2 result to chain file")
    //             .long_about(
    //                 r#"Convert minimap2 result to chain file

    // A paf file should be created with a command shown in below.

    // $ minimap2 -cx asm5 --cs NEW_FASTA ORIGINAL_FASTA > PAF_FILE.paf
    // "#,
    //             )
    //             .arg(
    //                 Arg::with_name("paf")
    //                     .index(1)
    //                     .takes_value(true)
    //                     .required(true)
    //                     .help("Input paf file"),
    //             )
    //             .arg(
    //                 Arg::with_name("chain")
    //                     .long("output")
    //                     .short("o")
    //                     .takes_value(true)
    //                     .required(true)
    //                     .help("Output chain file"),
    //             )
    //     }
    pub fn run(&self) -> anyhow::Result<()> {
        minimap2_to_chain_helper(&self.paf, &self.output)?;

        Ok(())
    }
}

// pub fn minimap2_to_chain(matches: &ArgMatches) {
//     match minimap2_to_chain_helper(
//         matches.value_of("paf").unwrap(),
//         matches.value_of("chain").unwrap(),
//     ) {
//         Ok(_) => (),
//         Err(e) => {
//             eprintln!("Error: {}", e);
//         }
//     }
// }

fn minimap2_to_chain_helper(paf_path: &str, chain_path: &str) -> anyhow::Result<()> {
    let mut paf_file = io::BufReader::new(open(paf_path).context("Cannot open paf file")?);
    let mut chain_file =
        create(chain_path, CompressionLevel::Default).context("Cannot create chain file")?;

    let mut count = 0;
    loop {
        count += 1;
        let mut line = String::new();
        if paf_file.read_line(&mut line)? == 0 {
            break; // finish
        }
        let elements: Vec<_> = line.trim().split('\t').collect();
        if elements.len() < 9 {
            return Err(anyhow::anyhow!("invalid format: line {}", count));
        }

        match elements[4] {
            "+" => {
                writeln!(
                    chain_file,
                    "chain\t4900\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    elements[0],
                    elements[1],
                    elements[2],
                    elements[3],
                    elements[5],
                    elements[6],
                    elements[4],
                    elements[7],
                    elements[8],
                    count
                )?;
            }
            "-" => {
                let seqlen = elements[6].parse::<u64>().expect("Cannot parse column 7");
                writeln!(
                    chain_file,
                    "chain\t4900\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    elements[0],
                    elements[1],
                    elements[2],
                    elements[3],
                    elements[5],
                    elements[6],
                    elements[4],
                    seqlen - elements[8].parse::<u64>().expect("Cannot parse column 9"),
                    seqlen - elements[7].parse::<u64>().expect("Cannot parse column 8"),
                    count
                )?;
            }
            _ => return Err(anyhow::anyhow!("invalid strand: line {}", count)),
        }

        static CIGER_MATCH: Lazy<Regex> = Lazy::new(|| Regex::new("(\\d+)([IMD])").unwrap());

        // parse CIGER
        let mut ciger = &elements
            .iter()
            .find(|x| x.starts_with("cg:Z:"))
            .expect("Cannot find CIGER")[5..];
        let mut ciger_list = Vec::new();

        while !ciger.is_empty() {
            let regexp_match = CIGER_MATCH
                .captures(ciger)
                .context(format!("Cannot parse ciger: line {}", count))?;
            if regexp_match.get(0).unwrap().start() != 0 {
                return Err(anyhow::anyhow!("Cannot parse: line {}", count));
            }
            //println!("{:?}", regexp_match);
            ciger_list.push((
                regexp_match
                    .get(1)
                    .unwrap()
                    .as_str()
                    .parse::<u32>()
                    .unwrap(),
                regexp_match
                    .get(2)
                    .unwrap()
                    .as_str()
                    .chars()
                    .next()
                    .unwrap(),
            ));
            ciger = &ciger[regexp_match.get(0).unwrap().end()..];
        }
        if elements[4] == "-" {
            ciger_list.reverse();
        }

        // Group CIGER
        let mut current_group = Vec::new();
        let mut chain_group = Vec::new();
        for one in ciger_list.iter() {
            if one.1 == 'M' && !current_group.is_empty() {
                chain_group.push(current_group.clone());
                current_group.clear()
            }
            current_group.push(one);
        }
        let last_group = current_group;

        // Write chain
        for one in chain_group {
            let matched: u32 = one.iter().filter(|x| x.1 == 'M').map(|x| x.0).sum();
            let insertion: u32 = one.iter().filter(|x| x.1 == 'I').map(|x| x.0).sum();
            let deletion: u32 = one.iter().filter(|x| x.1 == 'D').map(|x| x.0).sum();
            writeln!(chain_file, "{}\t{}\t{}", matched, insertion, deletion)?;
        }

        let matched: u32 = last_group.iter().filter(|x| x.1 == 'M').map(|x| x.0).sum();
        writeln!(chain_file, "{}", matched)?;
        writeln!(chain_file)?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use std::fs;

    #[test]
    fn test_minimap2chain() -> anyhow::Result<()> {
        fs::create_dir_all("../target/test-output/minimap2chain/")?;

        minimap2_to_chain_helper(
            "../liftover-rs/testfiles/genomes/chain/GRCh38-to-GRCh37.chr22.paf",
            "../target/test-output/minimap2chain/GRCh38-to-GRCh37.chain",
        )?;
        // TODO check result

        Ok(())
    }
}
