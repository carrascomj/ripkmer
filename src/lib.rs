use bio::io::fastq;
use bio_types::sequence::SequenceRead; // trait for len() of record
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;

pub struct Config {
    pub filename: String,
    pub k: usize,
    pub prefix: String,
}

impl Config {
    pub fn new(mut args: std::env::Args) -> Result<Config, &'static str> {
        args.next();

        let filename = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a filename!"),
        };

        let k = match args.next() {
            Some(arg) => arg.parse::<usize>().unwrap(),
            None => 16, // default
        };

        let prefix = match args.next() {
            Some(arg) => arg,
            None => String::from("ATCG"), // default
        };

        Ok(Config {
            filename,
            k,
            prefix,
        })
    }
}

/// Run the program
/// Read file -> Iterate over records -> Gather overlapping k-mers
///
/// # Errors
/// This function will return an error if it is incapable of reading the `filename`
pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let records = fastq::Reader::from_file(Path::new(&config.filename))?;
    let mut records = records.records();
    let all_kmers: HashMap<String, u32> = hash_kmer(&mut records, config.k, &config.prefix);
    // TODO: provisional
    println!("{:#?}", all_kmers);

    eprintln!("Running with {}-mers!", config.k);

    Ok(())
}

/// Apply kmerization over a whole file
fn hash_kmer(
    records: &mut fastq::Records<std::fs::File>,
    k: usize,
    prefix: &String,
) -> HashMap<String, u32> {
    let mut all_kmers: HashMap<String, u32> = HashMap::new();
    for record in records {
        let record = record.as_ref().unwrap();
        if record.check().is_ok() {
            kmerize(record, &mut all_kmers, k, prefix);
        }
    }
    all_kmers
}

/// Gather overlapping k-mers with a prefix
/// Usually just the k-mers starting with a prefix are used to reduce the size of the database.
fn kmerize(
    record: &fastq::Record,
    all_kmers: &mut HashMap<String, u32>,
    k: usize,
    prefix: &String,
) {
    let n = record.len();
    if n < k {
        // TODO: verify that this is the expected behaviour
        return;
    }

    let seq = record.seq();
    let prefix = prefix.as_bytes();
    for i in 0..(n - k + 1) {
        let kmer = &seq[i..(i + k)];
        if kmer.starts_with(prefix) {
            *all_kmers
                .entry(std::str::from_utf8(kmer).unwrap().to_string())
                .or_insert(0) += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_result() {
        let fq: &'static [u8] = b"@id description\nAATTAAGGAACC\n+\n!!!!!!!!!!!!\n";
        let records = fastq::Reader::new(fq)
            .records()
            .map(|record| record.unwrap());
        let mut all_kmers: HashMap<String, u32> = HashMap::new();

        let mut expected: HashMap<String, u32> = HashMap::new();
        expected.insert(String::from("AATT"), 1);
        expected.insert(String::from("AAGG"), 1);
        expected.insert(String::from("AACC"), 1);

        for record in records {
            kmerize(&record, &mut all_kmers, 4, &String::from("AA"))
        }
        assert_eq!(expected, all_kmers);
    }
    #[test]
    fn no_prefix() {
        let fq: &'static [u8] = b"@id description\nAATTA\n+\n!!!!!\n";
        let records = fastq::Reader::new(fq)
            .records()
            .map(|record| record.unwrap());

        let mut expected: HashMap<String, u32> = HashMap::new();
        expected.insert(String::from("AATT"), 1);
        expected.insert(String::from("ATTA"), 1);

        let mut all_kmers: HashMap<String, u32> = HashMap::new();
        for record in records {
            kmerize(&record, &mut all_kmers, 4, &String::from(""))
        }

        assert_eq!(expected, all_kmers);
    }
    #[test]
    fn n_less_than_k() {
        let fq: &'static [u8] = b"@id description\nAATTA\n+\n!!!!!\n";
        let records = fastq::Reader::new(fq)
            .records()
            .map(|record| record.unwrap());
        let mut all_kmers: HashMap<String, u32> = HashMap::new();
        for record in records {
            kmerize(&record, &mut all_kmers, 12, &String::from(""))
        }
        assert_eq!(HashMap::<String, u32>::new(), all_kmers);
    }
}
