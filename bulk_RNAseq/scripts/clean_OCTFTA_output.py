#!/usr/bin/python
# clean_OCTFTA_output.py
# Given the output of One Code To Find Them All (OCTFTA)
# remove those entries that are wrongly assembled.
# This script won't improve the OCTFTA assembly, but just
# clean it.
# How?
# - XXXX
# - YYYY
# - ZZZZ

# Input files:
# - TE GTF annotations with unique TE ID: 
#   /home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.gtf
# - danRer11.nonalt.fa.out.elem_sorted.csv
#   (Output of OCTFTA with all the TE annotations defragmented)
# - TE length file (2 columns, TE superfamily and TE length)
#??????

# Steps:
# - Read TE GTF into dictionary (save the original order)
# - Read TE length file and subset that all the TE 
#   superfamilies present on the TE GTF are in the file
# - Read OCTFTA file line by line and apply the following rules:
#   - If it's not an assembled TE, print the corresponding TE GTF line
#   - If it's an assembled TE:
#     - If 2 (or more) parts of the assemble map to the
#       same TE reference coordinates (X? of overlap), then
#       do NOT print the assembled TE but the corresponding TE GTF lines.
#     - If the length of the assembled TE is larger than
#       2x the length of the corresponding reference TE, 
#       do NOT print the assembled TE but the corresponding TE GTF lines.
#       (for LTRs, sum the LTR and Internal part)
#     - If everything is OK, print the assembled TE with 
#       a new unique TE ID, and the unique TE IDs of the
#       fragments assembled.


import argparse
import sys
import re
import os
from progress.bar import Bar
import time
from progress.bar import IncrementalBar

parser = argparse.ArgumentParser(description="Clean One Code To Find Them All output.")
parser.add_argument("-t", "--te_gtf", required = True, help="Path to TE gtf file.")
parser.add_argument("-i", "--input", required = True, help="Path to (OCTFTA) file.")
parser.add_argument("-l", "--te_length", required = True, help="Path to TE length file.")
parser.add_argument("-c", "--te_length_cutoff", required = False, default=1.5, type=float, help="TE reference sequence length fraction cutoff to discard large TEs. Value of 1 means that the cutoff is the length of the reference TE. Value of 2 means that the cutoff is 2x the length of the reference TE. Default: 1.5")
parser.add_argument("-o", "--output", required = True, help="Path to output file.")
args = vars(parser.parse_args())

TE_GTF_PATH = args['te_gtf']
OCTFTA_PATH = args['input']
TE_LEN_PATH = args['te_length']
TE_LEN_CUT  = args['te_length_cutoff']
OUTPUT_PATH = args['output']
#TE_GTF_PATH = "/home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.gtf"
#OCTFTA_PATH = "/home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/defragmented/danRer11.nonalt.fa.out.elem_sorted.csv"
#TE_LEN_PATH = "/home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/defragmented/tes.lengths"
#OUTPUT_PATH = "???"
#OCTFTA_PATH = "/home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/defragmented/danRer11.nonalt.fa.out.elem_sorted.head.csv"

### Functions
def parse_attribute_to_dict(attribute_string):
	attr_l = [x for x in [ i.lstrip() for i in attribute_string.split(";")] if x]
	attr_d = {d[0]:d[1].replace('"','') for d in [i.split(" ") for i in attr_l]}
	return(attr_d)

def unique_list(list1):
	# Function to get unique values from a list.
	# Copied from: https://www.geeksforgeeks.org/python-get-unique-values-list/#:~:text=Method%202%20%3A%20Using%20Set,a%20list%20to%20print%20it.
	# intilize a null list 
	unique_list = [] 
	# traverse for all elements 
	for x in list1: 
		# check if exists in unique_list or not 
		if x not in unique_list: 
			unique_list.append(x) 
	# Return 
	return unique_list

def get_TE_ref_seq_length(TEfrag_Ls, TElen_DICT):
	# Get the length of the TE reference sequence
	if len(set([i[9] for i in TEfrag_Ls])) == 1:
		#All elements in the TE frag lines are the same,
		#so only pick 1 TE from TElen_dict
		TE_ref_len = TElen_DICT[TEfrag_Ls[0][9]]
		return TE_ref_len
	else:
		#Different TEs in the frag lines, probably LTR-I-LTR.
		#Merge the length of each of the TEs.
		TE_ref_len = sum([TElen_dict[i[9]] for i in TEfrag_Ls])
		return TE_ref_len

def build_defrag_transcript_id(list_transcript_id):
	#defrag_transcript_id "L2-5_DRe_dup0/1/2"
	#defrag_transcript_id "EnSpm-9_DR_dup25/26/27"
	if False in ["_dup" in i for i in list_transcript_id]:
		# Some transcript_id do NOT have the "_dup" suffix
		TE_superfam = [i for i in list_transcript_id if "_dup" not in i][0]
		dup_array = "/".join([ re.sub(".+_dup","",i) for i in list_transcript_id if "_dup" in i ])
		defrag_transcript_id = "%s_dup0/%s" % (TE_superfam, dup_array)
	else:
		# All transcript_id have the "_dup" suffix
		TE_superfam = re.sub("_dup\d+","",list_transcript_id[0])
		dup_array = "/".join([ re.sub(".+_dup","",i) for i in list_transcript_id ])
		defrag_transcript_id = "%s_dup%s" % (TE_superfam, dup_array)
	#
	return defrag_transcript_id

def common_begining(str1,str2):
	# Check characters by order between str1 and str2
	# until they differ, and return the matching ones.
	match_chr = ""
	for i in range(0,min(len(str1),len(str2))):
		if str1[i] == str2[i]:
			match_chr += str1[i]
		else:
			return match_chr
			break
		#
	return match_chr

def build_defrag_transcript_id_v2(list_transcript_id):
	#First check if after removing "_dup\d+" from all TE ids,
	#if they are all the same. If that's the case, then it's easy.
	#If it's not the case, then for each different TE id 
	#(probably LTR and -I) join the _dup parts.
	if len(set( [re.sub("_dup\d+","",i) for i in list_transcript_id] )) == 1:
		#All TEs are the same, just with different _dupX
		if False in ["_dup" in i for i in list_transcript_id]:
			# Some transcript_id do NOT have the "_dup" suffix
			TE_superfam = [i for i in list_transcript_id if "_dup" not in i][0]
			dup_array = "/".join([ re.sub(".+_dup","",i) for i in list_transcript_id if "_dup" in i ])
			defrag_transcript_id = "%s_dup0/%s" % (TE_superfam, dup_array)
		else:
			# All transcript_id have the "_dup" suffix
			TE_superfam = re.sub("_dup\d+","",list_transcript_id[0])
			dup_array = "/".join([ re.sub(".+_dup","",i) for i in list_transcript_id ])
			defrag_transcript_id = "%s_dup%s" % (TE_superfam, dup_array)
	else:
		#Different TE IDs merged (probably LTR and I)
		diff_TEs = list(set( [re.sub("_dup\d+","",i) for i in list_transcript_id] ))
		common_str_TE = common_begining(diff_TEs[0],diff_TEs[1])
		# Note that if there's a 3rd diff_TE, I'm not checking it
		diff_TEs_defrag = []
		for dTE in diff_TEs:
			r = re.compile(dTE)
			sub_TEs = list(filter(r.match, list_transcript_id))
			if False in ["_dup" in i for i in sub_TEs]:
				# Some transcript_id do NOT have the "_dup" suffix
				TE_superfam = [i for i in sub_TEs if "_dup" not in i][0]
				dup_array = "/".join([ re.sub(".+_dup","",i) for i in sub_TEs if "_dup" in i ])
				sub_defrag_transcript_id = "%s_dup0/%s" % (TE_superfam, dup_array)
				sub_defrag_transcript_id = re.sub('/$','',sub_defrag_transcript_id)
				#print(sub_defrag_transcript_id)
			else:
				# All transcript_id have the "_dup" suffix
				TE_superfam = re.sub("_dup\d+","",sub_TEs[0])
				dup_array = "/".join([ re.sub(".+_dup","",i) for i in sub_TEs ])
				sub_defrag_transcript_id = "%s_dup%s" % (TE_superfam, dup_array)
				sub_defrag_transcript_id = re.sub('/$','',sub_defrag_transcript_id)
				#print(sub_defrag_transcript_id)
			#
			diff_TEs_defrag.append(sub_defrag_transcript_id)
		#
		defrag_transcript_id = common_str_TE+"_".join([T[len(common_str_TE):] for T in diff_TEs_defrag])
	#
	return defrag_transcript_id

def make_TEgtf_print_line(TEgtf_line):
	#i.e.: chr1	RepeatMasker	exon	469	611	494	+	.	gene_id "Harbinger-8N2_DR"; transcript_id "Harbinger-8N2_DR"; family_id "PIF-Harbinger"; class_id "DNA"; geneId "ENSDARG00000099104"; transcriptId "ENSDART00000158290"; annotation_fix "Distal Intergenic"; TE_GE_strand "OS"; annotation_2 "Intergenic.OS"; annotation_3 "Intergenic.OS"
	if 'defrag_transcript_id' not in TEgtf_line[8].keys():
		raise ValueError("Attribute defrag_transcript_id is missing")
	else:
		non_attribute = "\t".join(TEgtf_line[0:8])
		attribute_string = 'gene_id "%s"; transcript_id "%s"; defrag_transcript_id "%s"; family_id "%s"; class_id "%s"; geneId "%s"; transcriptId "%s"; annotation_fix "%s"; TE_GE_strand "%s"; annotation_2 "%s"; annotation_3 "%s"' % (TEgtf_line[8]['gene_id'], TEgtf_line[8]['transcript_id'], TEgtf_line[8]['defrag_transcript_id'], TEgtf_line[8]['family_id'], TEgtf_line[8]['class_id'], TEgtf_line[8]['geneId'], TEgtf_line[8]['transcriptId'], TEgtf_line[8]['annotation_fix'], TEgtf_line[8]['TE_GE_strand'], TEgtf_line[8]['annotation_2'], TEgtf_line[8]['annotation_3'])
		print_line = "%s\t%s" % (non_attribute, attribute_string)
		return print_line


### Run
# Read TE GTF
TEgtf_dict = {}
# Unnecessary # TEgtf_dict_chr = {}
TEgtf_dict_coord = {}
TEgtf_order = []
TEgtf_order_coord = []
print('Reading TE GTF file: %s' % (os.path.basename(TE_GTF_PATH)))
N_lines = sum(1 for line in open(TE_GTF_PATH))
bar = Bar('Reading', max=N_lines)
with open(TE_GTF_PATH, 'r') as FILE_IN:
	for line in FILE_IN:
		if line.startswith("#"):
			bar.next()
			continue
		#
		linesplit = line.strip().split("\t")
		attr_s = linesplit[8]
		attr_dict = parse_attribute_to_dict(attr_s)
		# Save it in the dict
		linesplit_sub = linesplit[:-1]
		linesplit_sub.append(attr_dict)
		TEgtf_dict[attr_dict['transcript_id']] = linesplit_sub
		# Save in coord dict
		#coord_string = "%s:%s-%s_%s_%s" % (linesplit[0],linesplit[3],linesplit[4],linesplit[6],attr_dict['gene_id'])
		coord_string = "%s:%s-%s_%s_%s_%s" % (linesplit[0],linesplit[3],linesplit[4],linesplit[6],linesplit[5],attr_dict['gene_id'])
		TEgtf_dict_coord[coord_string] = linesplit_sub
		# Save order in the list
		TEgtf_order.append(attr_dict['transcript_id'])
		# Save order in the list
		TEgtf_order_coord.append(coord_string)
		bar.next()

TEgtf_chroms = list(set([i.split(':')[0] for i in list(TEgtf_dict_coord.keys())]))
bar.finish()

# Read TE length file
TElen_dict = {}
print('Reading TE length file: %s' % (os.path.basename(TE_LEN_PATH)))
TE_superfam_uniq = unique_list([ value[8]['gene_id'] for key, value in TEgtf_dict.items() ])
N_lines = sum(1 for line in open(TE_LEN_PATH))
bar = Bar('Reading', max=N_lines)
with open(TE_LEN_PATH, 'r') as FILE_IN:
	for line in FILE_IN:
		#linesplit = line.strip().split("\t")# Some line have spaces instead of tab
		linesplit = line.strip().split()
		#Discard all TEs that are not in TE_superfam_uniq
		if linesplit[0] not in TE_superfam_uniq:
			bar.next()
			continue
		else:
			TElen_dict[linesplit[0]] = int(linesplit[1])
		bar.next()
bar.finish()


# Read OCTFTA file
print('Reading TE length file: %s' % (os.path.basename(OCTFTA_PATH)))
N_lines = sum(1 for line in open(OCTFTA_PATH))
bar = Bar('Reading', max=N_lines)
with open(OCTFTA_PATH, 'r') as FILE_IN, open(OUTPUT_PATH, 'w') as FILE_OUT:
	#Skip header
	header = FILE_IN.readline(); bar.next()
	#print(header.strip())
	for line in FILE_IN:
		if line == "\n":
			bar.next()
			continue
		if line == header:
			bar.next()
			continue
		#
		linesplit = line.strip().split("\t")
		linesplit[0] = linesplit[0].replace('###','')
		if linesplit[4] not in TEgtf_chroms:
			bar.next()
			continue
		#
		if "/" in linesplit[0]:
			#print("\"/\" in the first column.")
			#print(linesplit)
			TEdefrag_line = linesplit
			TEfrag_lines = []
			N_TEs = linesplit[0].count("/")
			for i in range(0, N_TEs+1):
				next_line = FILE_IN.readline(); bar.next()
				next_linesplit = next_line.strip().split("\t")
				TEfrag_lines.append(next_linesplit)
			#print(TEfrag_lines)
			#Check if the length of the assembled TE is larger than 1.5x the length of the reference TE sequence
			try:
				TE_ref_len = get_TE_ref_seq_length(TEfrag_lines, TElen_dict)
			except KeyError:
				#If it's not in the TE length dict, skip it.
				continue
			TE_defrag_len = int(TEdefrag_line[7])
			if TE_defrag_len > TE_ref_len*TE_LEN_CUT:
				#print("TE bigger than ref TE len x1.5")
				#Print frag TEs without joined defrag ID
				for TE_l in TEfrag_lines:
					try:
						#TEgtf_line = find_TE_in_TEgtf_dict_chrom(TE_l, TEgtf_dict_chr)
						coord_string = "%s:%s-%s_%s_%s_%s" % (TE_l[4],TE_l[5],TE_l[6],TE_l[8].replace("C","-"),TE_l[0],TE_l[9]) #coord_string = "%s:%s-%s_%s_%s_%s" % (linesplit[0],linesplit[3],linesplit[4],linesplit[6],linesplit[5],attr_dict['gene_id'])
						TEgtf_line = TEgtf_dict_coord[coord_string]
					except KeyError:
						print( "TE (%s) not found on GTF annotations (%s:%s-%s)" % (TE_l[9], TE_l[4],TE_l[5],TE_l[6]))
						continue
					#Print new TE GTF line
					defrag_transcript_id = TEgtf_line[8]['transcript_id']# Same as the original
					new_TEgtf_line = TEgtf_line
					new_TEgtf_line[8]['defrag_transcript_id'] = defrag_transcript_id
					print_TEgtf_line = make_TEgtf_print_line(new_TEgtf_line)
					FILE_OUT.write(print_TEgtf_line+"\n")# print(print_TEgtf_line)
			#Work in progress# Check if 2 (or more) parts of the TE fragments map the the same TE reference coordinates
			#Work in progress#   How?
			else:
				#print("TE smaller than ref TE len x1.5")
				# Print the frag TE with a joined defrag ID
				TEgtf_lines = []
				for TE_l in TEfrag_lines:
					try:
						#TEgtf_line = find_TE_in_TEgtf_dict_chrom(TE_l, TEgtf_dict_chr)
						coord_string = "%s:%s-%s_%s_%s_%s" % (TE_l[4],TE_l[5],TE_l[6],TE_l[8].replace("C","-"),TE_l[0],TE_l[9]) #coord_string = "%s:%s-%s_%s_%s_%s" % (linesplit[0],linesplit[3],linesplit[4],linesplit[6],linesplit[5],attr_dict['gene_id'])
						TEgtf_line = TEgtf_dict_coord[coord_string]
					except KeyError:
						print( "TE (%s) not found on GTF annotations (%s:%s-%s)" % (TE_l[9], TE_l[4],TE_l[5],TE_l[6]))
						continue
					#
					TEgtf_lines.append(TEgtf_line)
				#
				defrag_transcript_id = build_defrag_transcript_id_v2( [i[8]['transcript_id'] for i in TEgtf_lines] )
				for i in range(0,len(TEgtf_lines)):
					TEgtf_lines[i][8]['defrag_transcript_id'] = defrag_transcript_id
					print_TEgtf_line = make_TEgtf_print_line(TEgtf_lines[i])
					FILE_OUT.write(print_TEgtf_line+"\n")# print(print_TEgtf_line)
		#
		else:
			#print("\"/\" not in the first column.")
			#print(linesplit)
			# Find line in the TE gtf
			try:
				#TEgtf_line = find_TE_in_TEgtf_dict_chrom(linesplit, TEgtf_dict_chr)
				coord_string = "%s:%s-%s_%s_%s_%s" % (linesplit[4],linesplit[5],linesplit[6],linesplit[8].replace("C","-"),linesplit[0],linesplit[9]) #coord_string = "%s:%s-%s_%s_%s_%s" % (linesplit[0],linesplit[3],linesplit[4],linesplit[6],linesplit[5],attr_dict['gene_id'])
				TEgtf_line = TEgtf_dict_coord[coord_string]
			except KeyError:
				print( "TE (%s) not found on GTF annotations (%s:%s-%s)" % (linesplit[9], linesplit[4],linesplit[5],linesplit[6]))
				bar.next(); continue
			#Print new TE GTF line
			defrag_transcript_id = TEgtf_line[8]['transcript_id']# Same as the original
			new_TEgtf_line = TEgtf_line
			new_TEgtf_line[8]['defrag_transcript_id'] = defrag_transcript_id
			print_TEgtf_line = make_TEgtf_print_line(new_TEgtf_line)
			FILE_OUT.write(print_TEgtf_line+"\n")# print(print_TEgtf_line)
			bar.next()
bar.finish()

################################

# defrag_ids = ["L2-5_DRe","L2-5_DRe_dup1","L2-5_DRe_dup2"]
# defrag_ids = ["Helitron-2_DR_dup1","Helitron-2_DR_dup2","Helitron-2_DR_dup6"]
# defrag_ids = ["ERV1-N4-LTR_DR","ERV1-N4-I_DR","ERV1-N4-LTR_DR_dup1"]
# defrag_ids = ["BHIKHARI_LTR_dup3","BHIKHARI_I-int","BHIKHARI_I-int_dup1","BHIKHARI_LTR_dup4"]
# defrag_ids = ["Gypsy12-LTR_DR_dup13","Gypsy12-I_DR_dup26","Gypsy12-LTR_DR_dup14"]
