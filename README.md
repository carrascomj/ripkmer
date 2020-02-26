# ripkmer
[![Build Status](https://travis-ci.com/carrascomj/ripkmer.svg?branch=master)](https://travis-ci.com/carrascomj/ripkmer)  
There are two ways of viewing this:

- Some k-mer algorithms using [Rust-Bio](https://github.com/rust-bio/rust-bio/) [[1]](#koster2016).
- My first project in Rust just to get confident with it.

## Features
The first idea is to reproduce in Rust the [KmerFinder](https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/) [[2]](#kmerfinder2014)
(in Python, but also in [JavaScript](https://github.com/yosoyubik/kmerfinderjs-docker)).

* [x] K-mer count on FASTQ.
* [x] Filter by prefix.
* [ ] Make it work for FASTA and BED files.
* [x] Compare k-mer distribution of two inputs.
* [ ] Move towards a KMA implementation.

## CLI Example
For this example, the first two FASTQ files of 
[SRR396636](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR396636), corresponding
to reads from _Pseudomonas aeruginsa MPAO1/P1_, with 1909263 sequences of ~100 bp each, were downloaded.

Having **ripkmer** installed and in the `$PATH`:
```bash
ripkmer SRR396636.sra_1.fastq SRR396636.sra_2.fastq
```
where the `k` number and the `prefix` would be left as default, being equivalent
to:
```bash
ripkmer SRR396636.sra_1.fastq SRR396636.sra_2.fastq 16 ATGAC
```


The output is in tabular format and can be redirected to standard output (and should not take much more than 4s).

    (16-mers)	Unique	Redundant	Intersection_unique	Intersection
    SRR396636.sra_1.fastq	23196	97871	34.19%	58.81%
    SRR396636.sra_2.fastq	30698	89107	25.83%	64.59%

where
- **Unique** is the number of unique k-mers found in the file;
- **Redundant** is the number of total k-mers found (with repetitions);
- **Interesection_unique** is the number of common unique k-mers found in both files;
- and **Intersection** is the number of total common k-mers found.

## References

[<a name="koster2016">1</a>] Köster, J. (2016). Rust-Bio: a fast and safe bioinformatics library. Bioinformatics, 32(3), 444-446.  
[<a name="kmerfinder2014">2</a>] Benchmarking of Methods for Genomic Taxonomy. Larsen MV, Cosentino S, Lukjancenko O, Saputra D, Rasmussen S, Hasman H, Sicheritz-Pontén T, Aarestrup FM, Ussery DW, Lund O. J Clin Microbiol. 2014 Feb 26