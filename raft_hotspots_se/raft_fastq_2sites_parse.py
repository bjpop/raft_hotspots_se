#!/usr/bin/env python

'''
raft_fastq_2sites_parse.py

Initiated by DJP 28th Sep  2016
Takes fastq file from GAIIx single fragment read run generated from Tchurikov et al 2015 RAFT (rapid amplification of forum domains) procedure.
Produces a cleaned version so that the reads are trimmed of RAFT sequences and presented with the non-Sau3AI site at the 5-prime end of the read.
Depending on the insertion orientation during the library build, this can require reverse-complementing of the relevant portion of an originating read.
As a result, the double-strand break site appears at the 5-prime end of the read

DSB proximal
CCGAATTCTCCTTATACTGCAGGGG	double-strand break site	INSERT	Sau3AI site next 4 bases	GATCGTTTGCGGCCGCTTAAGCTTGGG

DSB distal and not necessarily contained within the read - depending on library insert size
CCCAAGCTTAAGCGGCCGCAAACGATC	Sau3AI site prev 4 bases	INSERT	double-strand break site	CCCCTGCAGTATAAGGAGAATTCGG

Version to avoid reading whole object into memory

Modified from raft_fastq_parse.py to require DSB and Sau3AI primers evidence. 
'''
from Bio import SeqIO as seqio
from Bio.Seq import Seq as seq
import argparse

parser = argparse.ArgumentParser(description='A tool to present single fragment fastq reads from RAFT trimmed with DSBs at 5-prime terminus')
parser.add_argument('--fastq',
                    metavar='FASTQ',
                    type=str,
                    required=True,
                    help='name of fastq file')
parser.add_argument('--dsb_flank',
                    metavar='DSB_FLANK',
                    type=str,
                    default='TGCAGGGG',
                    help='RAFT adapter sequence 5-prime proximal to the break site')
parser.add_argument('--re_flank',
                    metavar='RE_FLANK',
                    type=str,
                    default='GATCGTTT',
                    help='RAFT adapter sequence 3-prime proximal to the restriction endonuclease (e.g. Sau3AI) site')

def main():
    options = parser.parse_args()
    dsb_flank = options.dsb_flank
    re_flank = options.re_flank
    dsb_flank_rc = str(seq(dsb_flank).reverse_complement())
    re_flank_rc = str(seq(re_flank).reverse_complement())
    for read in seqio.parse(options.fastq, "fastq"):
	if dsb_flank in read.seq[:30]:		# checks for DSB proximal by end of RAFT primer in first 30 bases 
	    if (read.seq.count(dsb_flank) < 2):	# checks for lack of adapter concatamerisation
		if (dsb_flank_rc not in read.seq) and (re_flank_rc not in read.seq):	# checks for same primer at both ends and hetero-concatamerisation-type events
		    ins_start = str(read.seq).find(dsb_flank) + len(dsb_flank)	# grab the insert start index
		    if read.seq[30:].count(re_flank) == 1:	# checks for only one rc of Sau3AI primer after first 30 bases
			ins_end = str(read.seq).find(re_flank)	# grab the insert endnindex
			print '@'+str(read.description)
			print read[ins_start:ins_end].seq
			print '+'
			print ''.join([chr(x+33) for x in read[ins_start:ins_end].letter_annotations['phred_quality']])
	elif re_flank_rc in read.seq[:30]:		# checks for Sau3AI primer proximal
	    if read.seq.count(re_flank_rc) < 2:	# checks for lack of Sau3AI primer concatemerisation
		if (re_flank not in read.seq) and (dsb_flank not in read.seq):	# checks for Sau3AI primer both ends
		    ins_start = str(read.seq).find(re_flank_rc) + len(re_flank_rc)	# grab the insert start index
		    if read.seq[30:].count(dsb_flank_rc) == 1:	# checks for only one rc of non-Sau3AI-RAFT primer
			ins_end = str(read.seq).find(dsb_flank_rc)	# grab the insert endnindex
			print '@'+str(read.description)
			print read.seq[ins_start:ins_end].reverse_complement()
			print '+'
			print ''.join([chr(x+33) for x in read[ins_start:ins_end].letter_annotations['phred_quality'][::-1]])

if __name__ == '__main__':
    main()
