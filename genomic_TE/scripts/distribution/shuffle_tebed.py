#!/usr/bin/env python3

import random
import sys

with open(sys.argv[1]) as infile:
    data = [line.split() for line in infile]

shuffled_tes = [[i[3], i[4]] for i in data]

for i in range(1, 11):
    outfilename = f'../data/shuffled-beds/te_gene_shuffled_{i}.bed'
    random.shuffle(shuffled_tes)
    with open(outfilename, 'w') as outfile:
        for i, line in enumerate(data):
            newline = line[:3] + shuffled_tes[i] + line[5:]
            newline = '\t'.join(newline) + '\n'
            outfile.write(newline)

