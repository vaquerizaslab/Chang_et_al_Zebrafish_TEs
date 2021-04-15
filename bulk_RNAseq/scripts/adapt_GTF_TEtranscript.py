#!/usr/bin/python
# adapt_GTF_TEtranscript_fix_TEfam.py
# Adapt GTF file of TEs to the style of TEtranscripts.
# Include another input file with TEsubfamily/TEfamily/TEclass to fix 
# the annotation in the GTF to have always the same TE family/class.
# Also, if entry is not in the TE info file, then it won't be printed.
# Among the changes, set exon as the feature column instead of the TE subfamily name.

import argparse
from progress.bar import Bar
import sys

parser = argparse.ArgumentParser(description="Adapt GTF file of TEs to the style of TEtranscripts.")
parser.add_argument("-i", "--input", required = True, help="Input gtf file. Attribute fields TE_subfam, TE_fam and TE_class must be present.")
parser.add_argument("-c", "--classes", required = True, help="Input file with TE name/class/family info to correct the GTF if necessary. Should look like danRer11_classes_wredundant.subClass.txt")
parser.add_argument("-o", "--output", required = True, help="Path to the output GTF file")
args = vars(parser.parse_args())

INPUT_FILE = args['input']
TE_INFO_FILE = args['classes']
OUTPUT_FILE = args['output']

# Functions
def parse_attribute_to_dict(attribute_string):
	attr_l = [x for x in [ i.lstrip() for i in attribute_string.split(";")] if x]
	attr_d = {d[0]:d[1].replace('"','') for d in [i.split(" ") for i in attr_l]}
	return(attr_d)

# Read TE_INFO_FILE file first and build dictionary
#TE_INFO_FILE="/home/qrovira/Projects/ZF_GRCz11_TEs/data/19_11_19_TE_annotations_Jon/danRer11_classes_wredundant.subClass.txt"
TE_info_dict = {}
with open(TE_INFO_FILE, 'r') as FILE_IN:
	for line in FILE_IN:
		linesplit = line.strip().split("\t")
		TEname, TEclass, TEfamily = linesplit
		if TEname in TE_info_dict.keys():
			sys.exit("ERROR: --classes file must contain unique values in column 0 (TE subfamily)")
		else:
			#TE_info_dict[TEname] = [TEname, TEclass, TEfamily]
			TE_info_dict[TEname] = {'TE_class':TEclass, 'TE_fam':TEfamily}

# Run
TE_count = {}
N_lines_INPUT = sum(1 for line in open(INPUT_FILE))
bar = Bar('Adapting GTF', max=N_lines_INPUT)
with open(INPUT_FILE, 'r') as FILE_IN, open(OUTPUT_FILE, 'w') as FILE_OUT:
	for line in FILE_IN:
		linesplit = line.strip().split("\t")
		seqname, source, feature, start, end, score, strand, frame, attr_s = linesplit
		attr_dict = parse_attribute_to_dict(attr_s)
		# Check if TE is in TE_info dict
		if attr_dict['TE_subfam'] not in TE_info_dict.keys(): #Discard entry
			continue
		# Check if attr_dict matches TE info and change it if it doesn't
		if attr_dict['TE_class'] != TE_info_dict[attr_dict['TE_subfam']]['TE_class']:
			attr_dict['TE_class'] = TE_info_dict[attr_dict['TE_subfam']]['TE_class']
		if attr_dict['TE_fam'] != TE_info_dict[attr_dict['TE_subfam']]['TE_fam']:
			attr_dict['TE_fam'] = TE_info_dict[attr_dict['TE_subfam']]['TE_fam']
		# Count
		if attr_dict['TE_subfam'] not in TE_count.keys():
			TE_count[attr_dict['TE_subfam']] = 0
			TE_uniq_ID = "%s" % (attr_dict['TE_subfam'])
		else:
			TE_count[attr_dict['TE_subfam']] += 1
			TE_uniq_ID = "%s_dup%i" % (attr_dict['TE_subfam'], TE_count[attr_dict['TE_subfam']])
		# Print
		new_attr_s = 'gene_id "%s"; transcript_id "%s"; family_id "%s"; class_id "%s";' % (attr_dict['TE_subfam'], TE_uniq_ID, attr_dict['TE_fam'],attr_dict['TE_class'])
		new_line = [seqname, source, "exon", start, end, score, strand, frame, new_attr_s]
		FILE_OUT.write("\t".join(new_line)+"\n")
		bar.next()

bar.finish()
print("Done")
##########