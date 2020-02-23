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
    let records = records.records();
    let mut all_kmers: HashMap<String, u32> = HashMap::new();
    for record in records {
        let record = record.unwrap();
        if record.check().is_ok() {
            kmerize(record, &mut all_kmers, config.k, &config.prefix);
        }
    }
    // TODO: provisional
    println!("{:#?}", all_kmers);

    eprintln!("Running with {}-mers!", config.k);

    Ok(())
}

/// Gather overlapping k-mers with a prefix
/// Usually just the k-mers starting with a prefix are used to reduce the size of the database.
fn kmerize(record: fastq::Record, all_kmers: &mut HashMap<String, u32>, k: usize, prefix: &String) {
    let n = record.len();
    if n < k {
        // TODO: verify that this is the expected behaviour
        return;
    }
    let seq = std::str::from_utf8(&record.seq()).unwrap().to_string();
    for i in 0..(n - k + 1) {
        let kmer = String::from(&seq[i..(i + k)]);
        if kmer.starts_with(prefix) {
            *all_kmers.entry(kmer).or_insert(0) += 1;
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
            kmerize(record, &mut all_kmers, 4, &String::from("AA"))
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
            kmerize(record, &mut all_kmers, 4, &String::from(""))
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
            kmerize(record, &mut all_kmers, 12, &String::from(""))
        }
        assert_eq!(HashMap::<String, u32>::new(), all_kmers);
    }
}
