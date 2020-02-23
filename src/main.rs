use ripkmer::Config;
use std::env;
use std::process;

/// Find kmers starting with a given prefix in a file
/// # Usage
/// In the command-line:
///
/// ```bash
/// kmer filename.fastq k prefix
/// ```
///
/// Where `filename.fastq` is a path to a FASTQ file, `k` is the k-mer number and prefix is
/// as sequence to filter out k-mers that doesn't start with that sequence.
///
/// # Examples
/// In the command-line:
///
/// ```bash
/// kmer example.fastq 16 ATCG
/// ```
fn main() {
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    if let Err(e) = ripkmer::run(config) {
        eprintln!("Application error: {}", e);
        process::exit(2);
    }
}
