#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
from os.path import commonprefix
import sys

def build_ltrdict(filepath):
    """Assigns a common element name to each I-LTR pair."""
    ltrdict = {}
    with open(filepath) as infile:
        for line in infile:
            line = line.split()
            if len(line) == 2:
                ltrname = commonprefix(line)
                for suffix in ['-', '_DRe', '_DR', '_']:
                    ltrname = ltrname.rstrip(suffix)
            elif len(line) == 1:
                ltrname = line[0]
            for i in line:
                if ':' in i:
                    for j in i.split(':'):
                        ltrdict[j] = ltrname
                else:
                    ltrdict[i] = ltrname
    return ltrdict

def main():
    """Splits every TE family into fasta file containing defragged seqs."""
    
    defragdir = '../../data/repeatmasker-out/defragmented'
    seqdir = '../../data/seqs'
    minlen, minseqs, maxseqs = (int(i) for i in sys.argv[1:])

    ltrdict = build_ltrdict(f'{defragdir}/ltrdict.txt')

    c = 0
    splitrecords = defaultdict(list)
    for record in SeqIO.parse(f'{defragdir}/danRer11.nonalt.fa.out.tes.fasta',
                              'fasta'):
        element = record.id.split('|')[-1].split(':')[0]
        element = ltrdict.get(element, element)
        record.id = f'{element}|{c}'
        
        # Sensible min seq lenght is ~50.
        if len(record.seq) >= minlen:
            splitrecords[element].append(record)
            c += 1
    
    # Sort sequences by length and take t.
    # Ignore elements with less than 10 sequences.
    elementlist = []
    for element, records in splitrecords.items():
        if len(records) < minseqs:
            continue
        sorted_records = sorted(records, key=lambda x: len(x.seq), reverse=True)
        SeqIO.write(records[:maxseqs], f'{seqdir}/{element}.fa', 'fasta')
        elementlist.append(element)

    with open(f'{seqdir}/elementlist.txt', 'w') as outfile:
        for element in elementlist:
            outfile.write(f'{element}\n')
        
if __name__ == '__main__':
    # main()
    ltrdict = build_ltrdict('../../data/repeatmasker-out/defragmented/ltrdict.txt')
    for key, val in ltrdict.items():
        print(f'{key}\t{val}')
