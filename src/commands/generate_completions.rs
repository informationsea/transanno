use super::Command;
use crate::cli;
use clap::{crate_name, App, Arg, ArgMatches, Shell};
use std::fs;

pub struct GenerateCompletions;

impl Command for GenerateCompletions {
    fn command_name(&self) -> &'static str {
        "generate-completions"
    }
    fn config_subcommand(&self, app: App<'static, 'static>) -> App<'static, 'static> {
        app.about("Generate completion files").arg(
            Arg::with_name("output")
                .index(1)
                .help("Output directory")
                .required(true)
                .takes_value(true),
        )
    }
    fn run(&self, matches: &ArgMatches<'static>) -> Result<(), crate::LiftOverError> {
        let outdir = matches.value_of("output").unwrap();
        fs::create_dir_all(&outdir).expect("Cannot create output directory");
        let mut app = cli::build_cli();
        app.gen_completions(crate_name!(), Shell::Bash, &outdir);
        app.gen_completions(crate_name!(), Shell::Zsh, &outdir);
        app.gen_completions(crate_name!(), Shell::PowerShell, &outdir);
        app.gen_completions(crate_name!(), Shell::Fish, &outdir);
        Ok(())
    }
}
