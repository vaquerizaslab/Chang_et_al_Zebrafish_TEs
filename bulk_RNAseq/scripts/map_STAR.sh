#!/bin/bash
#map_STAR.sh
# Map reads with STAR

STAR_INDEX=$1
READS_PAIR_1=$2
READS_PAIR_2=$3
FILE_PREFIX=$4
THREADS=$5

# Run
echo -e "STAR --runThreadN $THREADS \\
     --genomeDir $STAR_INDEX \\
     --readFilesIn $READS_PAIR_1 $READS_PAIR_2 \\
     --readFilesCommand zcat \\
     --outSAMtype BAM Unsorted \\
     --chimSegmentMin 10 \\
     --winAnchorMultimapNmax 200 \\
     --outFilterMultimapNmax 100 \\
     --outFileNamePrefix $FILE_PREFIX"
STAR --runThreadN $THREADS \
     --genomeDir $STAR_INDEX \
     --readFilesIn $READS_PAIR_1 $READS_PAIR_2 \
     --readFilesCommand zcat \
     --outSAMtype BAM Unsorted \
     --chimSegmentMin 10 \
     --winAnchorMultimapNmax 200 \
     --outFilterMultimapNmax 100 \
     --outFileNamePrefix $FILE_PREFIX

##############