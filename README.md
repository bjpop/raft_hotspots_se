# Overview 

DNA breakage arises during a variety of biological processes, including transcription, replication and genome re-arrangements. 

A major advance in this direction comes from the work of Tchurikov et al, reporting that sites of human genome double-strand breaks (DSBs) occur frequently at sites in rDNA that are tightly linked with active transcription - the authors used a RAFT (rapid amplification of forum termini) protocol that selects for blunt-ended sites.

This suite of programs facilitates the analysis of the FASTQ outputs of the RAFT protocol, and allows the detection of double stranded breakpoint hotspots at a single nucleotide resolution.

This suite of tools was used in the work described by the paper
"Fine resolution mapping of double-strand break sites for human ribosomal DNA units" [Genomics Data 10 (2016) 19–21](http://www.sciencedirect.com/science/article/pii/S221359601630109X).

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bjpop/raft_hotspots_se/master/LICENSE)

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

We provide two Python programs (`raft_fastq_2sites_parse` and `raft_bed_2sites_parse`) which are intended to be used within a pipeline.
The pipeline requires bash, bwa (>= 0.7.5), samtools (>= 1.3.1), bedtools (>= 2.17.0), Python (2.7).

The pipeline supports the following arguments:

```
    pipeline.sh [-h] [-v] -r reference -f fastq -p prefix -b blacklist -s sites
```
  
* `-h`: (optional), print a help message and exit
* `-v`: (optional), print more information about progress to standard output
* `-f`: input FASTQ file
* `-p`: output filename prefix, to be used in creating new files
* `-b`: list of blacklisted genome sites to avoid, as a bed file
* `-s`: sites located within 5 base pairs of a Sau3AI consensus site (GATC), as a bed file

Example command line:

```
./pipeline.sh -v -r U13369.1_hg19.fa -f SRR944107.fastq -p SRRbla_hg19 -b hg19.blacklist.sort.bed -s hg19_GATC5.bed
```

We also provide an example SLURM script to illustrate running the example on a cluster.

The underlying Python scripts (`raft_fastq_2sites_parse` and `raft_bed_2sites_parse`) support the following command line arguments: 
```
% raft_fastq_2sites_parse -h
usage: raft_fastq_2sites_parse [-h] --fastq FASTQ [--dsb_flank DSB_FLANK]
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
% raft_bed_2sites_parse -h
usage: raft_bed_2sites_parse [-h] --bed BED

A tool to parse raft bamtobed output after presenting fastq reads with DSBs at
start. This yields DSB coordinates and counts.

optional arguments:
  -h, --help  show this help message and exit
  --bed BED   name of bed file
```

# Exit status values

Raft_hotspots_se returns the following exit status values:

* *0*: The program completed successfully.
* *1*: The program terminated with an error.
* *2*: There was an error with the command line arguments.


# Bugs

File at our [Issue Tracker](https://github.com/bjpop/raft_hotspots_se/issues)
