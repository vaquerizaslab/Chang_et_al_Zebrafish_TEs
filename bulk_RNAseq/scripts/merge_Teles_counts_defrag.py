#!/usr/bin/python
# merge_Teles_counts_defrag.py
# Given the join telescope results, merge the
# rows that contain TEs that have been joined
# in the defrag step, attribute defrag_transcript_id.

import argparse
import sys
import re
import os
from progress.bar import Bar
from more_itertools import unique_everseen

parser = argparse.ArgumentParser(description="Given the join telescope results, merge the rows that contain TEs that have been joined in the defrag step.")
parser.add_argument("--te_gtf", required = True, help="Path to TE GTF file. Must contain attribute field \"defrag_transcript_id\".")
parser.add_argument("--teles_join", required = True, help="Path to merged telescope count results.")
parser.add_argument("--output", required = True, help="Path to output file.")
parser.add_argument("--verbose", required = False, default = True, help="Print progress bar.")
args = vars(parser.parse_args())

TEGTF_PATH  = args['te_gtf']
TELES_PATH  = args['teles_join']
OUTPUT_PATH = args['output']
VERBOSE     = args['verbose']
#TEGTF_PATH  = "/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_09_13_TEdefrag_DET/TEdefrag_annot/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c1.gtf"
#TELES_PATH  = "/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_09_13_TEdefrag_DET/merge_Teles_counts/head.join_telesRep.final_count.tab"
#OUTPUT_PATH = "/home/qrovira/Projects/ZF_GRCz11_TEs/results/20_09_13_TEdefrag_DET/merge_Teles_counts/TEST.reads_FWD_tes_ALL.join_telesRep.final_count.dfrag.c1.tab"
#VERBOSE     = True

### Functions
def parse_attribute_to_dict(attribute_string):
	attr_l = [x for x in [ i.lstrip() for i in attribute_string.split(";")] if x]
	attr_d = {d[0]:d[1].replace('"','') for d in [i.split(" ") for i in attr_l]}
	return(attr_d)

### Run
# Read TE gtf file and make a dictionary of 
# transcript_id to defrag_transcript_id.
print( "Reading TE gtf file: %s" % (os.path.basename(TEGTF_PATH)) )
if VERBOSE==True:
	N_lines = sum(1 for line in open(TEGTF_PATH))
	bar = Bar('Reading', max=N_lines)
#TEGTF_PATH="/Users/quirze/cluster/Projects/ZF_GRCz11_TEs/results/20_09_13_TEdefrag_DET/TEdefrag_annot/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.head.gtf"
TE_2_TEdefrag_ID_dict = {}
#TEdefrag_2_TE_ID_dict = {}
#TEdefrag_ID_order = []
with open(TEGTF_PATH, 'r') as FILE_IN:
	for line in FILE_IN:
		if line.startswith("#"):
			if VERBOSE==True:
				bar.next()
			continue
		#
		linesplit = line.strip().split("\t")
		attr_s = linesplit[8]
		attr_dict = parse_attribute_to_dict(attr_s)
		#TEdefrag_ID_order.append(attr_dict['defrag_transcript_id'])
		#if attr_dict['defrag_transcript_id'] not in TEdefrag_ID_order:
		#	TEdefrag_ID_order.append(attr_dict['defrag_transcript_id'])
		#
		TE_2_TEdefrag_ID_dict[attr_dict['transcript_id']] = attr_dict['defrag_transcript_id']
		#try:
		#	TE_2_TEdefrag_ID_dict[attr_dict['transcript_id']].append(attr_dict['defrag_transcript_id'])
		#	#TEdefrag_2_TE_ID_dict[attr_dict['defrag_transcript_id']].append(attr_dict['transcript_id'])
		#except KeyError:
		#	TE_2_TEdefrag_ID_dict[attr_dict['transcript_id']] = [attr_dict['defrag_transcript_id']]
		#	#TEdefrag_2_TE_ID_dict[attr_dict['defrag_transcript_id']] = [attr_dict['transcript_id']]
		#
		if VERBOSE==True:
			bar.next()
if VERBOSE==True:
	bar.finish()

# Read the Telescope file and merge lines, when 
# necessary, and save into a dict.
print( "Reading Telescope counts file: %s" % (os.path.basename(TELES_PATH)) )
if VERBOSE==True:
	N_lines = sum(1 for line in open(TELES_PATH))
	bar = Bar('Reading', max=N_lines)
Teles_merge_dict = {}
TEdefrag_ID_order = []
with open(TELES_PATH, 'r') as FILE_IN:
	header = FILE_IN.readline()
	Teles_merge_dict['header'] = header.strip().split("\t")
	if VERBOSE==True:
		bar.next()
	#
	for line in FILE_IN:
		linesplit = line.strip().split("\t")
		TE_id = linesplit[0]
		if TE_id not in TE_2_TEdefrag_ID_dict.keys():
			if VERBOSE==True:
				bar.next()
			continue
		else:
			TEdefrag_id = TE_2_TEdefrag_ID_dict[TE_id]
			TEdefrag_ID_order.append(TEdefrag_id)
			try:
				Teles_merge_dict[TEdefrag_id] = [ x + y for x, y in zip(Teles_merge_dict[TE_id], [int(i.replace("NA","0")) for i in linesplit[1:]]) ]
			except KeyError:
				Teles_merge_dict[TEdefrag_id] = [int(i.replace("NA","0")) for i in linesplit[1:]]
			#
			if VERBOSE==True:
				bar.next()
if VERBOSE==True:
	bar.finish()

TEdefrag_ID_order_u = list(unique_everseen(TEdefrag_ID_order))

# Print output file
print( "Writing merged counts to: %s" % (os.path.basename(OUTPUT_PATH)) )
if VERBOSE==True:
	N_lines = len(TEdefrag_ID_order_u)
	bar = Bar('Writing', max=N_lines)
with open(OUTPUT_PATH, 'w') as FILE_OUT:
	#Header
	FILE_OUT.write('\t'.join(Teles_merge_dict['header'])+'\n')
	if VERBOSE==True:
		bar.next()
	#
	for TE in TEdefrag_ID_order_u:
		joined_counts = '\t'.join( [str(i) for i in Teles_merge_dict[TE]] )
		print_line = "%s\t%s" % (TE, joined_counts)
		FILE_OUT.write(print_line+"\n")
		#
		if VERBOSE==True:
			bar.next()
if VERBOSE==True:
	bar.finish()

print("** Done **")
##########