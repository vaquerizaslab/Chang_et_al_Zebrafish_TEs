#!/usr/bin/python
#join_telesRep_column.py
# Join the results of Telescope (after sorting them) into a single file.
# Given a list of files and a the column to use to merge, output the files stacked one next to each other.

import argparse
import sys

parser = argparse.ArgumentParser(description="Join the results of Telescope (after sorting them) into a single file. Given a list of files and a the column to use to merge, output the files stacked one next to each other.")
parser.add_argument("-f", "--files", required = True, help="List of files separated by comas \",\". i.e: File1.tab,File2.tab,File3.tab,...,FileN.tab ")
parser.add_argument("-c", "--column", required = True, type=int, help="Column to select when joining the files. This should be 0-based number.")
parser.add_argument("-n", "--names", required = True, help="List of names to use as header in the output file.")
parser.add_argument("-o", "--ouput", required = True, help="Output file name.")
args = vars(parser.parse_args())

FILES = args['files']
COLUMN  = args['column']
NAMES = args['names']
OUT_FILE = args['ouput']

# Read the files and save only the specified column
NAMES_list = NAMES.split(",")
FILES_list = FILES.split(",")
if len(NAMES_list) != len(FILES_list):
	sys.exit("ERROR: Length of files list and names list is different.")
#
ALL_columns = []
for F in FILES_list:
	print("Reading file: %s" % (F))
	with open(F, 'r') as FILE_IN:
		if ALL_columns == []:
			S_column = []
			TE_column = []
			for line in FILE_IN:
				if line.startswith('#'):
					continue
				linesplit = line.strip().split("\t")
				S_column.append(linesplit[COLUMN])
				TE_column.append(linesplit[0])
			ALL_columns.append(TE_column)
			ALL_columns.append(S_column)
		else:
			S_column = []
			for line in FILE_IN:
				if line.startswith('#'):
					continue
				linesplit = line.strip().split("\t")
				S_column.append(linesplit[COLUMN])
			ALL_columns.append(S_column)

# Now print the output
print("Writing file: %s" % (OUT_FILE))
with open(OUT_FILE, 'w') as FILE_OUT:
	FILE_OUT.write( "%s\t%s\n" % ("TEs", "\t".join(NAMES_list)) )
	for i in range(0,len(ALL_columns[0])):
		FILE_OUT.write( "%s\n" % ( "\t".join([j[i] for j in ALL_columns]) ) )

# Done
print("*** Done ***")
##################