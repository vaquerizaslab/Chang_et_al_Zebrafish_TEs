#!/bin/bash
#trim_PE_reads.sh
#Trim pair-end reads with bbduk

FASTQ1=$1
FASTQ2=$2
OUT_FASTQ1=$3
OUT_FASTQ2=$4
ADAPTERS="~/bin/bbmap/resources/adapters.fa"

echo "bbduk.sh in1=$FASTQ1 in2=$FASTQ2 out1=$OUT_FASTQ1 out2=$OUT_FASTQ2 ref=$ADAPTERS overwrite=t ktrim=r k=23 mink=11 hdist=1 tpe tbo"
bbduk.sh in1=$FASTQ1 in2=$FASTQ2 out1=$OUT_FASTQ1 out2=$OUT_FASTQ2 ref=$ADAPTERS overwrite=t ktrim=r k=23 mink=11 hdist=1 tpe tbo

############