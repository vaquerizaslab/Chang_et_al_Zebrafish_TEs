#!/bin/bash
#sambama_merge.sh
# First argument is the output BAM file
# Second argument is the number of CPUs to use
# The rest of the arguments are the BAM files to merge

ARGS=( $@ )
OUT_BAM=$1
THREADS=$2
TO_MERGE_BAMS=( ${ARGS[@]:2} )

echo "time sambamba merge -t $THREADS $OUT_BAM ${TO_MERGE_BAMS[*]}"
time sambamba merge -t $THREADS $OUT_BAM ${TO_MERGE_BAMS[*]}

#
echo "*** Done ***"
#############
