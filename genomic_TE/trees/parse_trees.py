#!/usr/bin/env python3

import os
import numpy as np
from scipy import stats
from ete3 import Tree
from Bio import SeqIO

def get_ancestors(tree):
    return list(set(leaf.up for leaf in tree.iter_leaves()))

def get_terminals(treefile):
    try:
        tetree = Tree(treefile)
    except:
        print(f'{treefile} could not be read')
        return
    ancestors = get_ancestors(tetree)
    erroneous = []
    for node in ancestors:
        pair = node.children
        dists = {pair[0].dist: pair[0], pair[1].dist: pair[1]}
        mind, maxd = min(dists.keys()), max(dists.keys())
        if (mind == 0 and maxd > 0.1):
            shorter = dists[mind]
            if shorter.name != '':
                erroneous.append(shorter.name)
        elif mind > 0 and maxd/mind > 100:
            shorter = dists[mind]
            if shorter.name != '':
                erroneous.append(shorter.name)
    prune_list = []
    for leaf in tetree.iter_leaves():
        if leaf.name not in erroneous:
            prune_list.append(leaf.name)
    tetree.prune(prune_list)
    return {leaf.name: leaf.dist for leaf in tetree.iter_leaves()}

def get_teclass(tename, seqdir):
    with open(f'{seqdir}/{tename}.fa') as infile:
        header = infile.readline()
    assert header.startswith(f'>{tename}')
    teclass = header.split()[1].split('|')[2].split('/')
    teclass, tefam = teclass[0], teclass[-1]
    return teclass, tefam

def get_fragment_lengths(tename, seqdir):
    fraglens = {}
    for record in SeqIO.parse(f'{seqdir}/{tename}.fa', 'fasta'):
        fraglens[record.id] = len(record.seq)
    return fraglens
    
def parse_trees(treedir, seqdir, outfilename):
    with open(outfilename, 'w') as outfile:
        outfile.write('element\ttename\ttefam\tteclass\tbranchlen\tfraglen\n')
        for filename in os.listdir(treedir):
            if not filename.endswith('.nwk'):
                continue
            terminals = get_terminals(f'{treedir}/{filename}')
            if not terminals:
                continue
            tename = filename.split('.')[0]
            teclass, tefam = get_teclass(tename, seqdir)
            fraglens = get_fragment_lengths(tename, seqdir)
            for el, blen in terminals.items():
                flen = fraglens[el]
                outfile.write(f'{el}\t{tename}\t{tefam}\t{teclass}\t{blen}\t{flen}\n')
        
def summarise_trees(treedir, seqdir, outfilename):
    with open(outfilename, 'w') as outfile:
        outfile.write('tename\ttefam\tteclass\tmeanlen\tmedlen\t25\t75\tvar\tnumbranches\n')
        for filename in os.listdir(treedir):
            if not filename.endswith('.nwk'):
                continue
            tename = filename.split('.')[0]
            teclass, tefam = get_teclass(tename, seqdir)
            
            terminals = get_terminals(f'{treedir}/{filename}')
            if not terminals:
                continue
            elif len(terminals) < 5:
                continue
            blens = list(terminals.values())
            desc = stats.describe(blens)
            bmean, bvar = desc.mean, desc.variance
            bmed = np.median(blens)
            b25, b75 = np.percentile(blens, 25), np.percentile(blens, 75)
            numbranches = len(blens)
            
            line = [tename, tefam, teclass, bmean, bmed, b25, b75, bvar,
                    numbranches]
            line = '\t'.join([str(i) for i in line])
            outfile.write(f'{line}\n')

def main():
    parse_trees('../../data/trees', 
                '../../data/seqs', 
                '../../data/trees/trees_all_repeats.txt')
    summarise_trees('../../data/trees', 
                    '../../data/seqs', 
                    '../../data/trees/trees_summary.txt')

if __name__ == '__main__':
    main()
    # get_terminals('../../data/trees/hAT-N126_DR.mafft.nwk')

