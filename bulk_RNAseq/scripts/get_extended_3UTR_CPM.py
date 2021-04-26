#!/usr/bin/python
# get_extended_3UTR_CPM.py

import argparse
import sys
import pyBigWig
from progress.bar import Bar

parser = argparse.ArgumentParser(description="Get mean and max value from given bigwig file for each region in given BED file.")
parser.add_argument("-i", "--input", required = True, help="Path to input BED file.")
parser.add_argument("-f", "--forward", required = True, help="Path to forward BigWig track.")
parser.add_argument("-r", "--reverse", required = True, help="Path to reverse BigWig track.")
parser.add_argument("-o", "--output", required = True, help="Path to output file name.")
args = vars(parser.parse_args())

BED_FILE = args['input']
BW_FWD_FILE  = args['forward']
BW_REV_FILE  = args['reverse']
OUT_FILE = args['output']

#
bw_fwd = pyBigWig.open(BW_FWD_FILE)
bw_rev = pyBigWig.open(BW_REV_FILE)
N_lines_INPUT = sum(1 for line in open(BED_FILE))
bar = Bar('Reading BED file and writing output', max=N_lines_INPUT)
with open(BED_FILE, 'r') as FILE_IN, open(OUT_FILE, 'w') as FILE_OUT:
	for line in FILE_IN:
		linesplit = line.strip().split("\t")
		chrom, start, end, Ensembl_gID, score, strand, gene_name, gene_biotype = linesplit
		if strand == "+":
			bw_mean = bw_fwd.stats(chrom, int(start), int(end), type="mean")[0]
			bw_max  = bw_fwd.stats(chrom, int(start), int(end), type="max")[0]
		elif strand == "-":
			bw_mean = bw_rev.stats(chrom, int(start), int(end), type="mean")[0]
			bw_max  = bw_rev.stats(chrom, int(start), int(end), type="max")[0]
		else:
			print("ERROR! Strand must be \"+\" or \"-\".")
		#
		FILE_OUT.write(line.strip() + '\t%.4f\t%.4f\n' % (bw_mean, bw_max) )
		#
		bar.next()

#########