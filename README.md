# Overview 

DNA breakage arises during a variety of biological processes, including transcription, replication and genome re-arrangements. 

A major advance in this direction comes from the work of Tchurikov et al, reporting that sites of human genome double-strand breaks (DSBs) occur frequently at sites in rDNA that are tightly linked with active transcription - the authors used a RAFT (rapid amplification of forum termini) protocol that selects for blunt-ended sites.

This suite of programs facilitates the analysis of the FASTQ outputs of the RAFT protocol, and allows the detection of double stranded breakpoint hotspots at a single nucleotide resolution.

This suite of tools was used in the work described by the paper
"Fine resolution mapping of double-strand break sites for human ribosomal DNA units"[Genomics Data](http://www.sciencedirect.com/science/article/pii/S221359601630109X).

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/raft_hotspots_se-paper/raft_hotspots_se/master/LICENSE)

# Installing

Raft_hotspots_se can be installed using `pip` in a variety of ways (`%` indicates the command line prompt):

1. Inside a virtual environment: 
```
% virtualenv raft_hotspots_se_dev
% source raft_hotspots_se_dev/bin/activate
% pip install -U /path/to/raft_hotspots_se
```
2. Into the global package database for all users:
```
% pip install -U /path/to/raft_hotspots_se
```
3. Into the user package database (for the current user only):
```
% pip install -U --user /path/to/raft_hotspots_se
```

# General behaviour


# Usage 

In the examples below, `%` indicates the command line prompt.

## Help message

```
% raft_fastq_parse -h
usage: raft_fastq_parse [-h] --fastq FASTQ [--dsb_flank DSB_FLANK]
                        [--re_flank RE_FLANK]

A tool to present single fragment fastq reads from RAFT trimmed with DSBs at
5-prime terminus

optional arguments:
  -h, --help            show this help message and exit
  --fastq FASTQ         name of fastq file
  --dsb_flank DSB_FLANK
                        RAFT adapter sequence 5-prime proximal to the break
                        site
  --re_flank RE_FLANK   RAFT adapter sequence 3-prime proximal to the
                        restriction endonuclease (e.g. Sau3AI) site
```

```
% raft_bed_parse -h
usage: raft_bed_parse [-h] --bed BED

A tool to parse raft bamtobed output after presenting fastq reads with DSBs at
start. This yields DSB coordinates and counts.

optional arguments:
  -h, --help  show this help message and exit
  --bed BED   name of bed file

```

# Exit status values

Raft_hotspots_se returns the following exit status values:

* *0*: The program completed successfully.
* *1*: File I/O error. This can occur if at least one of the input FASTA files cannot be opened for reading. This can occur because the file does not exist at the specified path, or raft_hotspots_se does not have permission to read from the file. 
* *2*: A command line error occurred. This can happen if the user specifies an incorrect command line argument. In this circumstance raft_hotspots_se will also print a usage message to the standard error device (stderr).
* *3*: Input FASTA file is invalid. This can occur if raft_hotspots_se can read an input file but the file format is invalid. 


# Bugs

File at our [Issue Tracker](https://github.com/bjpop/raft_hotspots_se/issues)
