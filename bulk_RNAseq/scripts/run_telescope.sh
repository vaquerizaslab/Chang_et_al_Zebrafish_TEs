#!/bin/bash
# run_telescope.sh

OUT_DIR=$1
SAMPLE=$2
BAM=$3
TE_ANNOT=$4

## Run
echo -e "telescope assign \\
  --theta_prior 200000 \\
  --max_iter 200 \\
  --attribute transcript_id \\
  --updated_sam \\
  --outdir $OUT_DIR \\
  --exp_tag telescope_${SAMPLE} \\
  $BAM $TE_ANNOT"
#
telescope assign \
  --theta_prior 200000 \
  --max_iter 200 \
  --attribute transcript_id \
  --updated_sam \
  --outdir $OUT_DIR \
  --exp_tag telescope_${SAMPLE} \
  $BAM $TE_ANNOT

echo "*** Done ***"
############