#!/usr/bin/python
# Convert RepeatMasker bed output (as received by Jon) to GTF format

import argparse
from progress.bar import Bar

parser = argparse.ArgumentParser(description="Convert RepeatMasker bed output (as received by Jon) to GTF format")
parser.add_argument("-i", "--input", required = True, help="Path to RepeatMasker output in BED format")
parser.add_argument("-o", "--output", required = True, help="Path to output file")
args = vars(parser.parse_args())

INPUT_FILE  = args['input']
OUTPUT_FILE = args['output']

# Functions
def bed_line_2_gtf(linesplit):
	# Get a list with the RM_bed line elements and return it as GTF.
	chrom, start, end, TEsubfam, SW_score, strand, perc_div, perc_del, perc_ins, pos_query_left, TEclassfam, pos_rep_begin, pos_rep_end, pos_rep_left, ID = linesplit[0:15]
	# Split TE class/fam
	if "/" in TEclassfam:
		TEclass, TEfam = TEclassfam.split("/")
	else:
		TEclass = TEclassfam
		TEfam = TEclassfam
	# Add +1 to start position
	start = str(int(start)+1)
	# Parse attributes
	attribute_string = 'TE_subfam "%s"; TE_fam "%s"; TE_class "%s"; SW_score "%s"; perc_div "%s"; perc_del "%s"; perc_ins "%s"; pos_query_left "%s"; pos_rep_begin "%s"; pos_rep_end "%s"; pos_rep_left "%s"; ID "%s"' % (TEsubfam, TEfam, TEclass, SW_score, perc_div, perc_del, perc_ins, pos_query_left, pos_rep_begin, pos_rep_end, pos_rep_left, ID)
	GTF_line = [chrom, "RepeatMasker", TEsubfam, start, end, SW_score, strand, ".", attribute_string]
	# Return
	return(GTF_line)

# Run
N_lines_INPUT = sum(1 for line in open(INPUT_FILE))

with open(INPUT_FILE, 'r') as FILE_IN, open(OUTPUT_FILE, 'w') as FILE_OUT:
	bar = Bar('Converting to GTF', max=N_lines_INPUT)
	for line in FILE_IN:
		linesplit = line.strip().split()
		line_gtf = bed_line_2_gtf(linesplit)
		FILE_OUT.write("\t".join(line_gtf)+"\n")
		bar.next()
	#
	bar.finish()

###########
