use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use quizx::json::{read_graph, write_graph};
use quizx::vec_graph::Graph;
use quizx_superopt::rewriter::CausalRewriter;
use quizx_superopt::superopt::{SuperOptOptions, SuperOptimizer};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// The input file
    #[arg(
        short,
        long,
        help = "The json-encoded pyzx graph file to use as input."
    )]
    input: PathBuf,
    /// The output file
    #[arg(
        short,
        long,
        help = "Output path for the optimized graph, encoded as pyzx json."
    )]
    output: PathBuf,
    /// The compiled rewriter.
    /// Defaults to `rewrites-2qb-lc.rwr` in the test files.
    #[arg(
        short,
        long,
        default_value = "../test_files/rewrites-2qb-lc.rwr",
        help = "The json-encoded pyzx graph file to use as input."
    )]
    rewriter: PathBuf,
    /// Timeout in seconds (default=no timeout)
    #[arg(
        short,
        long,
        value_name = "TIMEOUT",
        help = "Timeout in seconds (default=None)."
    )]
    timeout: Option<u64>,
    /// Maximum time in seconds to wait between circuit improvements (default=no timeout)
    #[arg(
        short = 'p',
        long,
        value_name = "PROGRESS_TIMEOUT",
        help = "Maximum time in seconds to wait between circuit improvements (default=None)."
    )]
    progress_timeout: Option<u64>,
    /// Max queue size.
    #[arg(
        short = 'q',
        long = "queue-size",
        default_value = "100",
        value_name = "QUEUE_SIZE",
        help = "The priority queue size. Defaults to 100."
    )]
    queue_size: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_time = Instant::now();

    // Parse command line arguments
    let args = Args::parse();

    println!("Loading input...");
    let graph: Graph = read_graph(&args.input).unwrap_or_else(|e| {
        panic!(
            "Could not read input graph {}.\n{e}",
            args.input.to_str().unwrap()
        )
    });

    println!("Loading rewrite set...");
    let f = File::open(&args.rewriter)?;
    let reader = BufReader::new(f);
    let rewriter: CausalRewriter<Graph> = rmp_serde::from_read(reader).unwrap_or_else(|e| {
        panic!(
            "Could not read compiled rewriter from {}.\n{e}",
            args.rewriter.to_str().unwrap()
        )
    });

    println!("Running the optimizer...");
    let optimizer: SuperOptimizer<CausalRewriter<Graph>> = SuperOptimizer::new(rewriter);
    let options = SuperOptOptions {
        timeout: args.timeout,
        progress_timeout: args.progress_timeout,
        queue_size: args.queue_size,
    };
    let result = optimizer.optimize(&graph, options);

    println!("Writing result...");
    write_graph(&result, true, &args.output)?;

    // Print the file size of output_file in megabytes
    let elapsed = start_time.elapsed();
    println!(
        "Done in {}.{:03} seconds",
        elapsed.as_secs(),
        elapsed.subsec_millis()
    );

    Ok(())
}
