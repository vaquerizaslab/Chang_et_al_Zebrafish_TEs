#!/usr/bin/env bash

# Same strand, upstream and downstream.
# Confirmed by manual inspection on genome browser
bedtools closest \
    -D 'b' \
    -s \
    -a ../data/repeatmasker-out/defragmented/danRer11.nonalt.fa.out.defrag.bed \
    -b ../data/genome/Danio_rerio.GRCz11.tss.bed |
    awk '{ OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $NF}' \
    > ../data/dist/GRCz11_defrag_te_closest_gene_ss_signed.bed

# Different strand, upstream and downstream.
bedtools closest \
    -D 'b' \
    -S \
    -a ../data/repeatmasker-out/defragmented/danRer11.nonalt.fa.out.defrag.bed \
    -b ../data/genome/Danio_rerio.GRCz11.tss.bed |
    awk '{ OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $NF}' \
    > ../data/dist/GRCz11_defrag_te_closest_gene_ds_signed.bed


