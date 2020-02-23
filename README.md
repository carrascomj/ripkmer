# ripkmer
There are two ways of viewing this:

- Some k-mer algorithms using [Rust-Bio](https://github.com/rust-bio/rust-bio/)[[1]](#amin2019).
- My first project in Rust just to get confident with it.

## Features
The first idea is to reproduce in Rust the [KmerFinder](https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/) [[2]](#amin2019)
(in Python, but also in [JavaScript](https://github.com/yosoyubik/kmerfinderjs-docker)).

* [x] K-mer count on FASTQ.
* [x] Filter by prefix.
* [ ] Make it work for FASTA and BED files.
* [ ] Compare k-mer distribution of two inputs.
* [ ] Move towards a KMA implementation.

## References

[<a name="koster2016">1</a>] Köster, J. (2016). Rust-Bio: a fast and safe bioinformatics library. Bioinformatics, 32(3), 444-446.  
[<a name="kmerfinder2014">2</a>] Benchmarking of Methods for Genomic Taxonomy. Larsen MV, Cosentino S, Lukjancenko O, Saputra D, Rasmussen S, Hasman H, Sicheritz-Pontén T, Aarestrup FM, Ussery DW, Lund O. J Clin Microbiol. 2014 Feb 26