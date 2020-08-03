use clap::{crate_authors, crate_version, App, Arg};

include!(concat!(env!("OUT_DIR"), "/git_version.rs"));

pub fn build_cli() -> App<'static, 'static> {
    App::new("transanno")
        .version(crate_version!())
        .long_version(concat!(crate_version!(), " git:", git_version!()))
        .author(crate_authors!())
        .about("Transfer annotation to other genome assemblies")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Verbose mode (-v: info, -vv: debug, -vvv: trace)"),
        )
        .subcommands(crate::commands::COMMANDS.iter().map(|x| x.cli()))
        .setting(clap::AppSettings::SubcommandRequiredElseHelp)
        .global_setting(clap::AppSettings::ColoredHelp)
        .global_setting(clap::AppSettings::ColorAuto)
}

pub fn validate_integer(text: String) -> Result<(), String> {
    for one in text.chars() {
        if !one.is_ascii_digit() {
            return Err(format!("{} is not number", text));
        }
    }
    Ok(())
}
