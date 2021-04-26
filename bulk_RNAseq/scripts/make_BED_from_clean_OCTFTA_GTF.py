#!/usr/bin/python
# make_BED_from_clean_OCTFTA_GTF.py
# Take the output of clean_OCTFTA_output.py
# (or a GTF file that has the attributes transcript_id 
# and defrag_transcript_id), and make a bed file of the
# defragmented TE "blocks"

import argparse
import sys
import re
import os
from progress.bar import Bar

parser = argparse.ArgumentParser(description="Make BED file from clean_OCTFTA_output.py output GTF.")
parser.add_argument("-g", "--gtf", required = True, help="Path to TE GTF file.")
parser.add_argument("-b", "--bed", required = True, help="Path to output BED file.")
args = vars(parser.parse_args())

GTF_PATH = args['gtf']
BED_PATH = args['bed']
#GTF_PATH = "/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_08_18_TEdefrag_overlap_DET/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.head.gtf"
#BED_PATH = "/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_08_18_TEdefrag_overlap_DET/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.head2.bed"

### Functions
def parse_attribute_to_dict(attribute_string):
	attr_l = [x for x in [ i.lstrip() for i in attribute_string.split(";")] if x]
	attr_d = {d[0]:d[1].replace('"','') for d in [i.split(" ") for i in attr_l]}
	return(attr_d)

def make_BED_line_from_GTF_line(GTF_list, GTF_attr_dict):
	# ==> danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.head.gtf <==
	# chr1	RepeatMasker	exon	1	397	1686	-	.	gene_id "Harbinger-N29_DR"; transcript_id "Harbinger-N29_DR"; defrag_transcript_id "Harbinger-N29_DR"; family_id "PIF-Harbinger"; class_id "DNA"; geneId "ENSDARG00000099104"; transcriptId "ENSDART00000158290"; annotation_fix "Distal"; TE_GE_strand "SS"; annotation_2 "Intergenic.SS"; annotation_3 "Intergenic.SS"
	# ==> danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.head.bed <==
	# chr1	0	397	Harbinger-N29_DR	1686	-	RepeatMasker	exon	.	gene_id "Harbinger-N29_DR"; transcript_id "Harbinger-N29_DR"; defrag_transcript_id "Harbinger-N29_DR"; family_id "PIF-Harbinger"; class_id "DNA"; geneId "ENSDARG00000099104"; transcriptId "ENSDART00000158290"; annotation_fix "Distal"; TE_GE_strand "SS"; annotation_2 "Intergenic.SS"; annotation_3 "Intergenic.SS"
	#
	seqname, source, feature, start, end, score, strand, frame = GTF_list[:8]
	BED_list = [seqname, int(start)-1, int(end), GTF_attr_dict['gene_id'], score, strand, source, feature, GTF_attr_dict['transcript_id'], GTF_attr_dict['defrag_transcript_id'], GTF_attr_dict['family_id'], GTF_attr_dict['class_id'], GTF_attr_dict['annotation_3']]
	BED_string = "\t".join([str(i) for i in BED_list]) 
	return(BED_string)

### Run
# Read TE GTF
N_lines = sum(1 for line in open(GTF_PATH))
bar = Bar('Reading', max=N_lines)
DEFRAG_lines_dict = {}
with open(GTF_PATH, 'r') as FILE_IN, open(BED_PATH, 'w') as FILE_OUT:
	for line in FILE_IN:
		if line.startswith("#"):
			bar.next()
			continue
		#
		linesplit = line.strip().split("\t")
		attr_s = linesplit[8]
		attr_dict = parse_attribute_to_dict(attr_s)
		#
		if attr_dict['transcript_id'] == attr_dict['defrag_transcript_id']:
			#Non defragmented TE. Print as BED
			BED_line = make_BED_line_from_GTF_line(linesplit[:8], attr_dict)
			FILE_OUT.write(BED_line+"\n")
		else:
			# Update the defrag start and end positions on the DEFRAG_lines_dict
			try:
				if int(linesplit[3]) < DEFRAG_lines_dict[attr_dict['defrag_transcript_id']][0][3]:
					DEFRAG_lines_dict[attr_dict['defrag_transcript_id']][0][3] = int(linesplit[3])
				if int(linesplit[4]) > DEFRAG_lines_dict[attr_dict['defrag_transcript_id']][0][4]:
					DEFRAG_lines_dict[attr_dict['defrag_transcript_id']][0][4] = int(linesplit[4])
			except KeyError:
				DEFRAG_lines_dict[attr_dict['defrag_transcript_id']] = [linesplit[:8], attr_dict]
				DEFRAG_lines_dict[attr_dict['defrag_transcript_id']][0][3] = int(DEFRAG_lines_dict[attr_dict['defrag_transcript_id']][0][3])
				DEFRAG_lines_dict[attr_dict['defrag_transcript_id']][0][4] = int(DEFRAG_lines_dict[attr_dict['defrag_transcript_id']][0][4])
			#
		#
		bar.next()
	#
	bar.finish()
	# Write the DEFRAG_lines_dict lines
	print("Writing defragmented lines")
	N_lines = len(DEFRAG_lines_dict)
	bar = Bar('Writing', max=N_lines)
	for key in list(DEFRAG_lines_dict.keys()):
		DEFRAG_linesplit = DEFRAG_lines_dict[key]
		BED_line = make_BED_line_from_GTF_line(DEFRAG_linesplit[0][:8], DEFRAG_linesplit[1])
		FILE_OUT.write(BED_line+"\n")
		bar.next()
	#
	bar.finish()

# Sort output BED file
print("Sorting output BED file")
print( "sort -k 1,1 -k2,2n %s > %s" % (BED_PATH, BED_PATH.replace(".bed", ".sorted.bed")) )
os.system( "sort -k 1,1 -k2,2n %s > %s" % (BED_PATH, BED_PATH.replace(".bed", ".sorted.bed")) )
print( "mv %s %s" % (BED_PATH.replace(".bed", ".sorted.bed"), BED_PATH) )
os.system( "mv %s %s" % (BED_PATH.replace(".bed", ".sorted.bed"), BED_PATH) )

##############