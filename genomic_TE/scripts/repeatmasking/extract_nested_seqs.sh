#!/usr/bin/env bash

ocdir='/usr/local/onecodetofindthemall'
rmaskdir='/Users/jonwells/Code/feschottelab/drerio-tes/data/repeatmasker-out'

# ltrdict.txt has been edited to resolve conflicting or ambiguous ltr names.
# tes.lengths are taken from repeatmasker database.
${ocdir}/one_code_to_find_them_all.pl \
    --rm ${rmaskdir}/danRer11.nonalt.fa.out \
    --ltr ${rmaskdir}/defragmented/ltrdict.txt \
    --length ${rmaskdir}/defragmented/tes.lengths \
    --fasta ${rmaskdir}/danRer11.nonalt.fa

# Combine output data into single files.
cat ${rmaskdir}/danRer11.nonalt.fa.out_*.ltr.csv > ${rmaskdir}/defragmented/danRer11.nonalt.fa.out.ltr.csv
cat ${rmaskdir}/danRer11.nonalt.fa.out_*.transposons.csv > ${rmaskdir}/defragmented/danRer11.nonalt.fa.out.transposons.csv
cat ${rmaskdir}/danRer11.nonalt.fa.out_*.elem_sorted.csv > ${rmaskdir}/defragmented/danRer11.nonalt.fa.out.elem_sorted.csv
cat ${rmaskdir}/danRer11.nonalt.fa.out_*.fasta > ${rmaskdir}/defragmented/danRer11.nonalt.fa.out.tes.fasta
${ocdir}/sum_copynumber.pl --dir ${rmaskdir} > ${rmaskdir}/defragmented/danRer11.nonalt.fa.out.copynumber.csv

# Clean up
rm ${rmaskdir}/danRer11.nonalt.fa.out_*.ltr.csv
rm ${rmaskdir}/danRer11.nonalt.fa.out_*.transposons.csv
rm ${rmaskdir}/danRer11.nonalt.fa.out_*.elem_sorted.csv
rm ${rmaskdir}/danRer11.nonalt.fa.out_*.copynumber.csv
rm ${rmaskdir}/danRer11.nonalt.fa.out_*.fasta
