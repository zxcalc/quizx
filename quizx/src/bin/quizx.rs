use clap::Parser;
use quizx::cli::Cli;

fn main() {
    let cli = Cli::parse();
    let result = cli.run();
    if let Err(e) = result {
        eprintln!("{e}");
        std::process::exit(1);
    }
}
