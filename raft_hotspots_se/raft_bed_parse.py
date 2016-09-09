#!/usr/bin/env python

'''
raft_bed_parse.py

Initiated by DJP 4th August  2016 - modified from raft_rdna_bed_parse.py
Takes bed file following raft_fastq_parse.py, bwa, sam to bam, sort, bamtobed including cigar.
Note that this processing has tried to clip off primer sequences and present the fastq reads so that the DSB site is at the 5' end of the read (start).
Note that during mapping, inserts could be oriented in either direction, so mapping can occur to + or - strands.
Further note that cigar relates to reference + strand. Minus strand mappings represent the reverse complement of a read mapped to + strand - the cigar reflects this. Plus strand mappings do not require transformation.

Example file lines:
reference	start	end	read	MAPQ	strand	cigar
gi|....|....	23	122	SRR944107.....	60	+	100M23S
chr1  386     421     SRR944107.26547584      47      +       6S35M
chr2  388     420     SRR944107.5704434       42      -       89S32M

We expect + strand mappings to have no clipping at the 'start'.
We expect - strand mappings to have to clipping at the 'end' - note again, this represents rc mapping.

Example command line:
./raft_bed_parse.py --bed <input.bed> > <output.bed>

3rd August 2016 - DJP modified for '-' strand to record base 0 instead of base 1.
3rd August 2016 - DJP modified to produce bed file format output for the first three columns followed by count

4th August 2016 - DJP modified to handle hg19-derived bed file with multiple chromosomes.
'''

import argparse

parser = argparse.ArgumentParser(description='A tool to parse raft bamtobed output after presenting fastq reads with DSBs at start.\
					      This yields DSB coordinates and counts.')
parser.add_argument('--bed',
                    metavar='BED',
                    type=str,
                    required=True,
                    help='name of bed file')

def get_cigar_l1(cigar):
    for character in cigar:
	if character.isdigit() != True:
	    return character

def main():
    options = parser.parse_args()
    dsb_site_dict = {}	# record sites by chr and, nested, co-ordinate (zero-based)
    chrom_set = set()	# record chr species encountered - a representation to be sorted prior to printout
    reader = open(options.bed)
    for read in reader:
	elements = read.split('\t')	# requires tab delimited format
	chrom, start, end = elements[0], int(elements[1]), int(elements[2])
	mq, strand, cigar = int(elements[4]), elements[5], elements[6]
	if (mq > 40) and ((end - start) > 25):	# check MAPQ > 40 and map length > 25
	    if strand == '+':
		if get_cigar_l1(cigar) != 'S':	# for + mapping check no clipping at 'start'
		    if chrom in chrom_set:
			if start in dsb_site_dict[chrom]:
			    dsb_site_dict[chrom][start] += 1
			else:
			    dsb_site_dict[chrom][start] = 1
		    else:
			dsb_site_dict[chrom] = {}
			dsb_site_dict[chrom][start] = 1
		    chrom_set.add(chrom)
		else:
		    continue
	    elif strand == '-':
		if cigar[-1] != 'S':	# for - mapping check no clipping at 'end'
		    end_base_zero = (end - 1)	# convert to base 0
		    if chrom in chrom_set:
			if end_base_zero in dsb_site_dict[chrom]:
			    dsb_site_dict[chrom][end_base_zero] += 1
			else:
			    dsb_site_dict[chrom][end_base_zero] = 1
		    else:
			dsb_site_dict[chrom] = {}
			dsb_site_dict[chrom][end_base_zero] = 1
		    chrom_set.add(chrom)
		else:
		    continue
    sorted_chroms = sorted(chrom_set)
    for contig in sorted_chroms:
	sorted_coords = sorted(dsb_site_dict[contig].keys())
	for coord in sorted_coords:
	    print '\t'.join(str(item) for item in [contig, coord, (coord + 1), dsb_site_dict[contig][coord]])

if __name__ == '__main__':
    main()
