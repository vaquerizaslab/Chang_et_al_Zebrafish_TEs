#!/usr/bin/env python3

import sys
import re
import pandas as pd
from collections import defaultdict
import csv

def get_familydict():
    with open('../../data/repeatmasker-out/danRer11_classes.txt') as infile:
        header = infile.readline()
        familydict = {l.split()[0]: l.split()[1:] for l in infile}
    return familydict

def reformat_attribute(x, familydict, dupcount):
    gene_id = re.search(r'"Motif:(.+)"', x).group(1)
    transcript_id = gene_id + f'_dup{dupcount[gene_id]}'
    family_id = familydict[gene_id][0]
    class_id = familydict[gene_id][1]
    dupcount[gene_id] += 1
    return f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; family_id "{family_id}"; class_id "{class_id}";'

def main():
    df = pd.read_csv(sys.argv[1], 
                     sep='\t', 
                     comment='#',
                     header=None,
                     names=['seqname', 'source', 'feature', 'start', 'end',
                            'score', 'strand', 'frame', 'attribute'])
    familydict = get_familydict()
    dupcount = defaultdict(lambda: 1)
    df['attribute'] = df['attribute'].apply(lambda x: reformat_attribute(x, familydict, dupcount))
    df['score'] = '.'
    df.to_csv(sys.argv[2], 
              index=False, 
              header=None, 
              sep='\t', 
              quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
    
