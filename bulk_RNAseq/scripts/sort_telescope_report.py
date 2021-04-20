#!/usr/bin/python
#sort_telescope_report.py
# Use GTF TE annotation to sort the output report of telescope

import argparse

parser = argparse.ArgumentParser(description="Use GTF TE annotation to sort the output report of telescope")
parser.add_argument("-g", "--gtf", required = True, help="GTF TE annotation")
parser.add_argument("-r", "--report", required = True, help="Telescope report file")
parser.add_argument("-o", "--ouput", required = True, help="Output file name")
args = vars(parser.parse_args())

GTF_FILE = args['gtf']
REPORT_FILE  = args['report']
OUT_FILE = args['ouput']

# Read the telescope report file
print("Reading the telescope report file")
Telescope_report_dict = {}
Telescope_report_stats = ""
with open(REPORT_FILE,'r') as FILE_IN:
	header = False
	for line in FILE_IN:
		if line.startswith('#'):
			Telescope_report_stats = line.strip()
		#
		linesplit = line.strip().split("\t")
		if header == False:
			header = linesplit
		else:
			transcript, transcript_length, final_count, final_conf, final_prop, init_aligned, unique_count, init_best, init_best_random, init_best_avg, init_prop = linesplit
			Telescope_report_dict[transcript] = linesplit

# Read GTF file. If TE is in telescope report then plot the results of telescope report.
# Otherwise, print all 0s or NA.
print("Reading the TE annotation file")
with open(GTF_FILE,'r') as FILE_IN, open(OUT_FILE,'w') as FILE_OUT :
	FILE_OUT.write(Telescope_report_stats+'\n')
	for line in FILE_IN:
		if line.startswith('#'):
			continue
		#
		linesplit = line.strip().split("\t")
		seqname, source, feature, start, end, score, strand, frame, attributes = linesplit
		attributes_dict = {i.replace('"','').split(" ")[0]: i.replace('"','').split(" ")[1] for i in attributes.split('; ')}
		TE_ID = attributes_dict['transcript_id']
		#
		try:
			FILE_OUT.write("\t".join(Telescope_report_dict[TE_ID])+'\n')
		except KeyError:
			FILE_OUT.write("%s\t%s\n" % (TE_ID, "\t".join(["NA"]*10)))

# Done
print("*** Done ***")
####################