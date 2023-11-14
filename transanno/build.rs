use std::env;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::process;
use std::str;

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("git_version.rs");
    let mut f = File::create(&dest_path).unwrap();

    let git_revision = match process::Command::new("git")
        .arg("rev-parse")
        .arg("HEAD")
        .output()
    {
        Ok(output) => {
            if output.stdout.is_empty() {
                "unknown".to_string()
            } else {
                str::from_utf8(&output.stdout).unwrap().trim().to_string()
            }
        }
        Err(_err) => "unknown".to_string(),
    };

    let git_status = match process::Command::new("git")
        .arg("status")
        .arg("--short")
        .output()
    {
        Ok(output) => {
            if output.stdout.is_empty() {
                ""
            } else {
                "-modified"
            }
        }
        Err(_err) => "-unknown",
    };

    //writeln!(f, "const GIT_REVISION : &str = \"{}\";", git_revision);
    //writeln!(f, "const GIT_VERSION : &str = \"{}{}\";", &git_revision[0..7], git_status);

    writeln!(
        f,
        r#"
macro_rules! git_version {{
    () => ("{}{}")
}}
"#,
        &git_revision[..7],
        git_status
    )
    .unwrap();
}
