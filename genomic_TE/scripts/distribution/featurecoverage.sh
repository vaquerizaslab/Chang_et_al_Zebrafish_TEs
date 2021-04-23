#!/usr/bin/env bash

genomedir="../data/genome"
lauedir="../data/Laue2019/"
repeatdir="../data/repeatmasker-out"
circosdir="../data/circosfig"

tewindow=2000000
genewindow=5000000
geneoverlap=3000000

# Make windows
echo "making windows"
bedtools makewindows \
    -g ${genomedir}"/Danio_rerio.GRCz11.karyotype.txt" \
    -w $tewindow > ${genomedir}"/tewindow.bed"

bedtools makewindows \
    -g ${genomedir}"/Danio_rerio.GRCz11.karyotype.txt" \
    -w $genewindow \
    -s $geneoverlap > ${genomedir}"/genewindow.bed"

# Make gene and chromatin tracks
echo "calculating gene coverage"
bedtools coverage \
    -a ${genomedir}"/tewindow.bed" \
    -b ${lauedir}"/GSM3096184_Sorted_SizedMG_6_2_H3K9me3-116183126.R1_001_trimmedMarkedDups.bam.DupsR.bedgraph" |
    awk 'BEGIN {OFS="\t"} { print $1, $2, $3, $7}' |
    cut -c 4- > ${circosdir}"/H3K9me3.density.txt"

bedtools coverage \
    -a ${genomedir}"/tewindow.bed" \
    -b ${genomedir}"/Danio_rerio.GRCz11.genes.bed" |
    awk 'BEGIN {OFS="\t"} { print $1, $2, $3, $7}' |
    cut -c 4- > ${circosdir}"/gene.density.txt"

# # Make TE tracks
# echo "calculating TE coverage"
bedtools coverage \
    -a ${genomedir}"/tewindow.bed" \
    -b ${repeatdir}"/danRer11.nonalt.fa.SINE.young1.bed" |
    awk 'BEGIN {OFS="\t"} { print $1, $2, $3, $7}' |
    cut -c 4- > ${circosdir}"/SINE.density.young1.txt"

bedtools coverage \
    -a ${genomedir}"/tewindow.bed" \
    -b ${repeatdir}"/danRer11.nonalt.fa.RC.young1.bed" |
    awk 'BEGIN {OFS="\t"} { print $1, $2, $3, $7}' |
    cut -c 4- > ${circosdir}"/RC.density.young1.txt"

bedtools coverage \
    -a ${genomedir}"/tewindow.bed" \
    -b ${repeatdir}"/danRer11.nonalt.fa.LINE.young1.bed" |
    awk 'BEGIN {OFS="\t"} { print $1, $2, $3, $7}' |
    cut -c 4- > ${circosdir}"/LINE.density.young1.txt"

bedtools coverage \
    -a ${genomedir}"/tewindow.bed" \
    -b "${repeatdir}/danRer11.nonalt.fa.LTR.young1.bed" |
    awk 'BEGIN {OFS="\t"} { print $1, $2, $3, $7}' |
    cut -c 4- > ${circosdir}"/LTR.density.young1.txt"

bedtools coverage \
    -a ${genomedir}"/tewindow.bed" \
    -b ${repeatdir}"/danRer11.nonalt.fa.DNA.young1.bed" |
    awk 'BEGIN {OFS="\t"} { print $1, $2, $3, $7}' |
    cut -c 4- > ${circosdir}"/DNA.density.young1.txt"

bedtools coverage \
    -a ${genomedir}"/tewindow.bed" \
    -b ${repeatdir}"/danRer11.nonalt.fa.BRSATI.bed" |
    awk 'BEGIN {OFS="\t"} { print $1, $2, $3, $7}' |
    cut -c 4- > ${circosdir}"/BRSATI.density.txt"
# Build circos figure
circos -conf ${circosdir}"/circos.conf"

