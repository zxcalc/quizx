use std::{fs, io::Write, time::Instant};

use clap::Parser;
use quizx::vec_graph::Graph;
use quizx_superopt::{rewrite_sets::RewriteSet, rewriter::CausalRewriter};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// The input file
    #[arg(short, long)]
    input: String,
    /// The output file
    #[arg(short, long)]
    output: String,
}

fn main() {
    let start_time = Instant::now();

    // Parse command line arguments
    let Args { input, output } = Args::parse();
    let inf = fs::File::open(&input).unwrap();

    println!("Loading rewrite rules from {:?}...", input);
    let rewrite_rules: Vec<RewriteSet<Graph>> = serde_json::from_reader(inf).unwrap();

    println!("Compiling {} rewrite rule sets...", rewrite_rules.len());
    let rewriter = CausalRewriter::from_rewrite_rules(rewrite_rules);

    println!("Writing rewriter to {:?}...", output);
    let mut outf = fs::File::create(&output).unwrap();
    outf.write(&rmp_serde::to_vec(&rewriter).unwrap()).unwrap();

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
}
