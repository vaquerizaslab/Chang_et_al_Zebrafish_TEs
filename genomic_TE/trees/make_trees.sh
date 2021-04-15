#!/usr/bin/env bash

repmaskdir='../../data/repeatmasker-out/defragmented'
seqdir='../../data/seqs'
treedir='../../data/trees'

################################################################################
## Extract all element families from defragmented RepeatMasker output
################################################################################

minlen=100
minseqs=10
maxseqs=1250

echo "min sequence length: ${minlen}"
echo "min num sequences: ${minseqs}"
echo "max num sequences: ${maxseqs}"
echo "splitting seqs..."

./split_defragmented_seqs.py $minlen $minseqs $maxseqs

echo "split seqs."

################################################################################
## Align sequences
################################################################################

echo "Aligning sequences..."

for te in `cat "${seqdir}/elementlist.txt"`; do
    /programs/mafft/bin/mafft \
        --thread 20 \
        --quiet \
        --auto "${seqdir}/${te}.fa" > "${seqdir}/${te}.mafft.fa"
done

echo "Aligned many a sequence."

################################################################################
## Remove highly gapped regions
################################################################################

echo "Cleaning gaps..."

for te in `cat "${seqdir}/elementlist.txt"`; do
   /local/workdir/jnw72/Software/trimal-1.4.1/source/trimal \
       -in "${seqdir}/${te}.mafft.fa" \
       -out "${seqdir}/${te}.tmp.mafft.fa" \
       -gt 0.01
   mv "${seqdir}/${te}.tmp.mafft.fa" "${seqdir}/${te}.mafft.fa"
done

echo "Cleaned gaps."

################################################################################
## Build trees
################################################################################

echo "Building trees..."

for te in `cat "${seqdir}/elementlist.txt"`; do
    numseqs=`grep '>' "${seqdir}/${te}.mafft.fa" -c`
    if [ $numseqs -gt 50 ]; then
        /programs/FastTree-2.1.10/FastTree -nt -gtr -pseudo -gamma -nopr\
            "${seqdir}/${te}.mafft.fa" > "../../data/trees/${te}.mafft.nwk"
    else
        /programs/FastTree-2.1.10/FastTree -nt -gtr -pseudo -nopr\
            "${seqdir}/${te}.mafft.fa" > "../../data/trees/${te}.mafft.nwk"
    fi
done

echo "Built many a tree."

