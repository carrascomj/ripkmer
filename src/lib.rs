use bio::io::fastq;
use bio_types::sequence::SequenceRead; // trait for len() of record
use std::cmp;
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;

pub struct Config {
    pub filename: String,
    pub filedb: String,
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
        let filedb = match args.next() {
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
            filedb,
            k,
            prefix,
        })
    }
}

/// Get the number of kmers (unique and total)
/// Useful for displaying and debugging.
struct Kstats {
    kunique: usize,
    kredundant: u32,
}

impl Kstats {
    fn new(kmers: &HashMap<String, u32>) -> Kstats {
        let mut kunique = 0;
        let mut kredundant = 0;
        for (_, count) in kmers {
            kunique += 1;
            kredundant += count;
        }
        Kstats {
            kunique,
            kredundant,
        }
    }
}

impl std::fmt::Display for Kstats {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "unique k-mers: {}, redundant k-mers: {}\n",
            self.kunique, self.kredundant
        )
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
    let kmers_target: HashMap<String, u32> = hash_kmer(&mut records, config.k, &config.prefix);

    let records = fastq::Reader::from_file(Path::new(&config.filedb))?;
    let mut records = records.records();

    let kmers_ref: HashMap<String, u32> = hash_kmer(&mut records, config.k, &config.prefix);
    let kstat_target = Kstats::new(&kmers_target);
    let kstat_ref = Kstats::new(&kmers_ref);
    let match_unique = intersect_keys(&kmers_target, &kmers_ref);
    let match_redundant = intersect_counters(&kmers_target, &kmers_ref);

    // TODO: this should be moved to Kstats somehow
    println!(
        "({}-mers)\tUnique\tRedundant\tIntersection_unique\tIntersection",
        config.k
    );
    println!(
        "{}\t{}\t{}\t{:.2}%\t{:.2}%",
        config.filename,
        kstat_target.kunique,
        kstat_target.kredundant,
        100f64 * match_unique as f64 / kstat_target.kunique as f64,
        100f64 * match_redundant as f64 / kstat_target.kredundant as f64 
    );
    println!(
        "{}\t{}\t{}\t{:.2}%\t{:.2}%",
        config.filedb,
        kstat_ref.kunique,
        kstat_ref.kredundant,
        100f64 * match_unique as f64 / kstat_ref.kunique as f64,
        100f64 * match_redundant as f64/ kstat_ref.kredundant as f64
    );

    Ok(())
}

/// Count intersection of keys in two HashMaps
fn intersect_keys<K: Eq + std::hash::Hash, V1, V2>(
    left: &HashMap<K, V1>,
    right: &HashMap<K, V2>,
) -> usize {
    left.keys()
        .filter(|k| right.contains_key(k))
        .count()
}

/// Count intersection of appeareances in two Counters (HashMaps)
fn intersect_counters<K: Eq + std::hash::Hash>(
    left: &HashMap<K, u32>,
    right: &HashMap<K, u32>,
) -> u32 {
    left.keys()
        .filter(|k| right.contains_key(k))
        .map(|k| cmp::min(left.get(k), right.get(k)).unwrap())
        .sum()
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
    #[test]
    fn kmer_stats() {
        let fq: &'static [u8] = b"@id description\nAATTAAGGAACC\n+\n!!!!!!!!!!!!\n";
        let mut all_kmers: HashMap<String, u32> = HashMap::new();
        let records = fastq::Reader::new(fq)
            .records()
            .map(|record| record.unwrap());
        for record in records {
            kmerize(&record, &mut all_kmers, 2, &String::from(""))
        }
        let calculated_kstats = Kstats::new(&all_kmers);
        
        assert_eq!((9, 11), (calculated_kstats.kunique, calculated_kstats.kredundant));
    }
    #[test]
    fn counter_comp() {
        let mut counter1: HashMap<String, u32> = HashMap::new();
        counter1.insert(String::from("a"), 23);
        counter1.insert(String::from("b"), 2);
        counter1.insert(String::from("c"), 15);
        let mut counter2: HashMap<String, u32> = HashMap::new();
        counter2.insert(String::from("a"), 5);
        counter2.insert(String::from("b"), 7);
        counter2.insert(String::from("c"), 3);
        assert_eq!(3, intersect_keys(&counter1, &counter2));
        assert_eq!(10, intersect_counters(&counter1, &counter2));
    }
}
