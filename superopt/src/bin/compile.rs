use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::{fs, io::Write, time::Instant};

use clap::Parser;
use quizx::vec_graph::Graph;
use quizx_superopt::{rewrite_sets::RewriteSet, rewriter::CausalRewriter};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// The input file
    #[arg(short = 'i', long)]
    input: PathBuf,
    /// The output file
    #[arg(short = 'o', long)]
    output: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_time = Instant::now();

    // Parse command line arguments
    let Args { input, output } = Args::parse();

    println!("Loading rewrite rules from {:?}...", input);
    let f = fs::File::open(&input)?;
    let reader = BufReader::new(f);
    let rewrite_rules: Vec<RewriteSet<Graph>> = serde_json::from_reader(reader)?;

    println!("Compiling {} rewrite rule sets...", rewrite_rules.len());
    let rewriter = CausalRewriter::from_rewrite_rules(rewrite_rules);

    println!("Writing rewriter to {:?}...", output);
    let f = fs::File::create(&output)?;
    let mut writer = BufWriter::new(f);
    writer.write(&rmp_serde::to_vec(&rewriter)?)?;

    // Print the file size of output_file in megabytes
    if let Ok(metadata) = fs::metadata(&output) {
        let file_size = metadata.len() as f64 / (1024.0 * 1024.0);
        println!("File size: {:.2} MB", file_size);
    }
    let elapsed = start_time.elapsed();
    println!(
        "Done in {}.{:03} seconds",
        elapsed.as_secs(),
        elapsed.subsec_millis()
    );

    Ok(())
}
