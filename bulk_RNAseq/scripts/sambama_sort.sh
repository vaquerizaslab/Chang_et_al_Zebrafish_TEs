#!/bin/bash
#sambama_sort.sh

BAM_FILE=$1
THREADS=$2

# sambamba generates the index (.bai) automatically
echo "sambamba sort -t ${THREADS} -o ${BAM_FILE/.bam/_sorted.bam} ${BAM_FILE}"
sambamba sort -t ${THREADS} -o ${BAM_FILE/.bam/_sorted.bam} ${BAM_FILE}
echo "mv ${BAM_FILE/.bam/_sorted.bam} ${BAM_FILE}"
mv ${BAM_FILE/.bam/_sorted.bam} ${BAM_FILE}
echo "mv ${BAM_FILE/.bam/_sorted.bam.bai} ${BAM_FILE}.bai"
mv ${BAM_FILE/.bam/_sorted.bam.bai} ${BAM_FILE}.bai

#############