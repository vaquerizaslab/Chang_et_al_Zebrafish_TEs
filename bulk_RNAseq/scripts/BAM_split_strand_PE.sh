#!/bin/bash
#BAM_split_strand_PE.sh

# Here we used the script to split BAM files by strand from Joseph Chao-Chung Kuo in his blog: 
# https://josephcckuo.wordpress.com/2016/11/18/splitting-reads-in-a-bam-file-by-strands/

# Get the bam file from the command line
BAM=$1
DIR=$2
THREADS=$3

FILE=$(basename $BAM)
NAME=${FILE%.*}
BAMF1=${DIR}/${NAME}_fwd1.bam
BAMF2=${DIR}/${NAME}_fwd2.bam
BAMF=${DIR}/${NAME}_fwd.bam
BAMR1=${DIR}/${NAME}_rev1.bam
BAMR2=${DIR}/${NAME}_rev2.bam
BAMR=${DIR}/${NAME}_rev.bam

# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse strand
#
# 0x1 - paired
# 0x2 - properly paired
# 0x20 - partner on reverse strand
# 0x40 - read one
# FLAGs 0x1 + 0x2 + 0x20 + 0x40 = 0x63 = 99 in decimal
echo "samtools view -bh -f 99 -@ $THREADS $BAM > $BAMF1"
samtools view -bh -f 99 -@ $THREADS $BAM > $BAMF1
echo "samtools index -@ $THREADS $BAMF1"
samtools index -@ $THREADS $BAMF1
# 0x1 - paired
# 0x2 - properly paired
# 0x10 - on reverse strand
# 0x80 - read two
# FLAGs 0x1 + 0x2 + 0x10 + 0x80 = 0x93 = 147 in decimal
echo "samtools view -bh -f 147 -@ $THREADS $BAM > $BAMF2"
samtools view -bh -f 147 -@ $THREADS $BAM > $BAMF2
echo "samtools index -@ $THREADS $BAMF2"
samtools index -@ $THREADS $BAMF2

#
# Combine alignments that originate on the forward strand.
#
echo "samtools merge -f -@ $THREADS $BAMF $BAMF1 $BAMF2"
samtools merge -f -@ $THREADS $BAMF $BAMF1 $BAMF2
echo "samtools index -@ $THREADS $BAMF"
samtools index -@ $THREADS $BAMF

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#

# 0x1 - paired
# 0x2 - properly paired
# 0x10 - reverse strand
# 0x40 - read one
# FLAGs 0x1 + 0x2 + 0x10 + 0x40 = 0x53 = 83 in decimal
echo "samtools view -bh -f 83 -@ $THREADS $BAM > $BAMR1"
samtools view -bh -f 83 -@ $THREADS $BAM > $BAMR1
echo "samtools index -@ $THREADS $BAMR1"
samtools index -@ $THREADS $BAMR1
# 0x1 - paired
# 0x2 - properly paired
# 0x30 - partner on reverse strand
# 0x80 - read two
# FLAGs 0x1 + 0x2 + 0x20 + 0x80 = 0xA3 = 163 in decimal
echo "samtools view -bh -f 163 -@ $THREADS $BAM > $BAMR2"
samtools view -bh -f 163 -@ $THREADS $BAM > $BAMR2
echo "samtools index -@ $THREADS $BAMR2"
samtools index -@ $THREADS $BAMR2

#
# Combine alignments that originate on the reverse strand.
#
echo "samtools merge -f -@ $THREADS $BAMR $BAMR1 $BAMR2"
samtools merge -f -@ $THREADS $BAMR $BAMR1 $BAMR2
echo "samtools index -@ $THREADS $BAMR"
samtools index -@ $THREADS $BAMR

echo "*** Done ***"
####################