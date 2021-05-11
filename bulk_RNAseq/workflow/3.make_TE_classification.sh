#!/bin/bash
# 3.make_TE_classification.sh

### Get de-fragmented annotations
# To reconstructs fragmented repeats and full-length LTR elements from TE annotations we used the
# tool onecodetofindthemall.pl from Bailly-Bechet et al., 2014. DOI: https://doi.org/10.1186/1759-8753-5-13
# The output of this tool should be named danRer11.nonalt.fa.out.elem_sorted.csv
# and placed in the ./data directory.
# A 2 column file with the length of the consensus sequence for each TE is also necessary: ./data/tes.lengths

# We then clean the output of onecodetofindthemall.pl to discard large TEs that have been joined
# based on the reference TE length.

TE_GTF="./data/danRer11.TEtrans_uID.gtf"
time python ./scripts/clean_OCTFTA_output.py \
  --te_gtf ./data/danRer11.TEtrans_uID.gtf \
  --input ./data/danRer11.nonalt.fa.out.elem_sorted.csv \
  --te_length ./data/tes.lengths \
  --te_length_cutoff 2 \
  --output ./data/danRer11.TEtrans_uID.dfrag.gtf
# Transform to BED
time python ./scripts/make_BED_from_clean_OCTFTA_GTF.py \
  --gtf ./data/danRer11.TEtrans_uID.dfrag.gtf \
  --bed ./data/danRer11.TEtrans_uID.dfrag.bed



### Merge BAMs
# In order to have more sequencing depth, replicates were merged before transcriptome assembly.
BAM_DIR="./data/BAMs"
BAM_DIR_MERGE="./data/BAMs/merged"
THREADS="3"
mkdir $BAM_DIR_MERGE
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/1_cell.merge.bam           $THREADS $BAM_DIR/1_cell_rep_1.bam           $BAM_DIR/1_cell_rep_2.bam           $BAM_DIR/1_cell_rep_3.bam           $BAM_DIR/1_cell_rep_4.bam           $BAM_DIR/1_cell_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/2_cell.merge.bam           $THREADS $BAM_DIR/2_cell_rep_1.bam           $BAM_DIR/2_cell_rep_2.bam           $BAM_DIR/2_cell_rep_3.bam           $BAM_DIR/2_cell_rep_4.bam           $BAM_DIR/2_cell_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/128_cell.merge.bam         $THREADS $BAM_DIR/128_cell_rep_1.bam         $BAM_DIR/128_cell_rep_2.bam         $BAM_DIR/128_cell_rep_3.bam         $BAM_DIR/128_cell_rep_4.bam         $BAM_DIR/128_cell_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/1k_cell.merge.bam          $THREADS $BAM_DIR/1k_cell_rep_1.bam          $BAM_DIR/1k_cell_rep_2.bam          $BAM_DIR/1k_cell_rep_3.bam          $BAM_DIR/1k_cell_rep_4.bam          $BAM_DIR/1k_cell_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/dome.merge.bam             $THREADS $BAM_DIR/dome_rep_1.bam             $BAM_DIR/dome_rep_2.bam             $BAM_DIR/dome_rep_3.bam             $BAM_DIR/dome_rep_4.bam             $BAM_DIR/dome_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/50pc_epiboly.merge.bam     $THREADS $BAM_DIR/50pc_epiboly_rep_1.bam     $BAM_DIR/50pc_epiboly_rep_2.bam     $BAM_DIR/50pc_epiboly_rep_3.bam     $BAM_DIR/50pc_epiboly_rep_4.bam     $BAM_DIR/50pc_epiboly_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/shield.merge.bam           $THREADS $BAM_DIR/shield_rep_1.bam           $BAM_DIR/shield_rep_2.bam           $BAM_DIR/shield_rep_3.bam           $BAM_DIR/shield_rep_4.bam           $BAM_DIR/shield_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/75pc_epiboly.merge.bam     $THREADS $BAM_DIR/75pc_epiboly_rep_1.bam     $BAM_DIR/75pc_epiboly_rep_2.bam     $BAM_DIR/75pc_epiboly_rep_3.bam     $BAM_DIR/75pc_epiboly_rep_4.bam     $BAM_DIR/75pc_epiboly_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/1_4_somites.merge.bam      $THREADS $BAM_DIR/1_4_somites_rep_1.bam      $BAM_DIR/1_4_somites_rep_2.bam      $BAM_DIR/1_4_somites_rep_3.bam      $BAM_DIR/1_4_somites_rep_4.bam      $BAM_DIR/1_4_somites_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/14_19_somites.merge.bam    $THREADS $BAM_DIR/14_19_somites_rep_1.bam    $BAM_DIR/14_19_somites_rep_2.bam    $BAM_DIR/14_19_somites_rep_3.bam    $BAM_DIR/14_19_somites_rep_4.bam    $BAM_DIR/14_19_somites_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/20_25_somites.merge.bam    $THREADS $BAM_DIR/20_25_somites_rep_1.bam    $BAM_DIR/20_25_somites_rep_2.bam    $BAM_DIR/20_25_somites_rep_3.bam    $BAM_DIR/20_25_somites_rep_4.bam    $BAM_DIR/20_25_somites_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/prim_5.merge.bam           $THREADS $BAM_DIR/prim_5_rep_1.bam           $BAM_DIR/prim_5_rep_2.bam           $BAM_DIR/prim_5_rep_3.bam           $BAM_DIR/prim_5_rep_4.bam           $BAM_DIR/prim_5_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/prim_15.merge.bam          $THREADS $BAM_DIR/prim_15_rep_1.bam          $BAM_DIR/prim_15_rep_2.bam          $BAM_DIR/prim_15_rep_3.bam          $BAM_DIR/prim_15_rep_4.bam          $BAM_DIR/prim_15_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/prim_25.merge.bam          $THREADS $BAM_DIR/prim_25_rep_1.bam          $BAM_DIR/prim_25_rep_2.bam          $BAM_DIR/prim_25_rep_3.bam          $BAM_DIR/prim_25_rep_4.bam          $BAM_DIR/prim_25_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/long_pec.merge.bam         $THREADS $BAM_DIR/long_pec_rep_1.bam         $BAM_DIR/long_pec_rep_2.bam         $BAM_DIR/long_pec_rep_3.bam         $BAM_DIR/long_pec_rep_4.bam         $BAM_DIR/long_pec_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/protruding_mouth.merge.bam $THREADS $BAM_DIR/protruding_mouth_rep_1.bam $BAM_DIR/protruding_mouth_rep_2.bam $BAM_DIR/protruding_mouth_rep_3.bam $BAM_DIR/protruding_mouth_rep_4.bam $BAM_DIR/protruding_mouth_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/day_4.merge.bam            $THREADS $BAM_DIR/day_4_rep_1.bam            $BAM_DIR/day_4_rep_2.bam            $BAM_DIR/day_4_rep_3.bam            $BAM_DIR/day_4_rep_4.bam            $BAM_DIR/day_4_rep_5.bam
bash ./scripts/sambama_merge.sh $BAM_DIR_MERGE/day_5.merge.bam            $THREADS $BAM_DIR/day_5_rep_1.bam            $BAM_DIR/day_5_rep_2.bam            $BAM_DIR/day_5_rep_3.bam            $BAM_DIR/day_5_rep_4.bam            $BAM_DIR/day_5_rep_5.bam

# Make transcriptome coverage tracks that will be used later
# First, split merged BAM files by strand get coverage tracks for forward and reverse strand
BAM_DIR_MERGE="./data/BAMs/merged"
BAM_DIR_MERGE_SPLIT="./data/BAMs/merged/strand_split"
mkdir $BAM_DIR_MERGE_SPLIT
BAM_FILES=( $BAM_DIR_MERGE/*.bam )
THREADS=6
for BAM in ${BAM_FILES[*]}; do
    bash ./scripts/BAM_split_strand_PE.sh $BAM $BAM_DIR_MERGE_SPLIT $THREADS
done
# Now, make tracks with bamCoverage
BAM_DIR_MERGE_ALL="./data/BAMs/merged"
BAM_DIR_MERGE_STRAND="./data/BAMs/merged/strand_split"
STAGES=( $(lh $BAM_DIR_MERGE_ALL | grep "bam" | rev | cut -d' ' -f1 | rev | cut -d'.' -f1 | sort -u) )
BW_DIR="./data/BWs"
mkdir $BW_DIR
for S in ${STAGES[*]}; do
    BAM_ALL="${BAM_DIR_MERGE_ALL}/${S}.merge.bam"
    BAM_FWD="${BAM_DIR_MERGE_STRAND}/${S}.merge_fwd.bam"
    BAM_REV="${BAM_DIR_MERGE_STRAND}/${S}.merge_rev.bam"
    #
    bamCoverage -p 4 --binSize 1 --verbose --minMappingQuality 10 --normalizeUsing CPM --samFlagInclude 64 -b $BAM_ALL -of bigwig -o $BW_DIR/${S}.CPM.bw
    bamCoverage -p 2 --binSize 1 --verbose --minMappingQuality 10 --normalizeUsing CPM --samFlagInclude 64 -b $BAM_FWD -of bigwig -o $BW_DIR/${S}.fwd.CPM.bw
    bamCoverage -p 2 --binSize 1 --verbose --minMappingQuality 10 --normalizeUsing CPM --samFlagInclude 64 -b $BAM_REV -of bigwig -o $BW_DIR/${S}.rev.CPM.bw
done



### Run stringtie
STR_DIR="./data/stringtie"
BAM_DIR_MERGE="./data/BAMs/merged"
THREADS="5"
mkdir $STR_DIR
stringtie $BAM_DIR_MERGE/1_cell.merge.bam           -o $STR_DIR/1_cell.stringtie.gtf           -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/2_cell.merge.bam           -o $STR_DIR/2_cell.stringtie.gtf           -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/128_cell.merge.bam         -o $STR_DIR/128_cell.stringtie.gtf         -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/1k_cell.merge.bam          -o $STR_DIR/1k_cell.stringtie.gtf          -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/dome.merge.bam             -o $STR_DIR/dome.stringtie.gtf             -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/50pc_epiboly.merge.bam     -o $STR_DIR/50pc_epiboly.stringtie.gtf     -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/shield.merge.bam           -o $STR_DIR/shield.stringtie.gtf           -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/75pc_epiboly.merge.bam     -o $STR_DIR/75pc_epiboly.stringtie.gtf     -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/1_4_somites.merge.bam      -o $STR_DIR/1_4_somites.stringtie.gtf      -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/14_19_somites.merge.bam    -o $STR_DIR/14_19_somites.stringtie.gtf    -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/20_25_somites.merge.bam    -o $STR_DIR/20_25_somites.stringtie.gtf    -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/prim_5.merge.bam           -o $STR_DIR/prim_5.stringtie.gtf           -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/prim_15.merge.bam          -o $STR_DIR/prim_15.stringtie.gtf          -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/prim_25.merge.bam          -o $STR_DIR/prim_25.stringtie.gtf          -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/long_pec.merge.bam         -o $STR_DIR/long_pec.stringtie.gtf         -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/protruding_mouth.merge.bam -o $STR_DIR/protruding_mouth.stringtie.gtf -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/day_4.merge.bam            -o $STR_DIR/day_4.stringtie.gtf            -p $THREADS --rf -t -c 1.5 -f 0.05
stringtie $BAM_DIR_MERGE/day_5.merge.bam            -o $STR_DIR/day_5.stringtie.gtf            -p $THREADS --rf -t -c 1.5 -f 0.05

# Get extended 3' UTRs for each sample
STR_DIR="./data/stringtie"
EXT_3UTR_DIR="./data/stringtie/ext_3UTR"
GE_GTF="./data/Danio_rerio.GRCz11.98.chr.ERCC.addChr.gtf.gz"
mkdir $EXT_3UTR_DIR
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/1_cell.stringtie.gtf           -r $GE_GTF -o $EXT_3UTR_DIR/1_cell.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/2_cell.stringtie.gtf           -r $GE_GTF -o $EXT_3UTR_DIR/2_cell.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/128_cell.stringtie.gtf         -r $GE_GTF -o $EXT_3UTR_DIR/128_cell.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/1k_cell.stringtie.gtf          -r $GE_GTF -o $EXT_3UTR_DIR/1k_cell.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/dome.stringtie.gtf             -r $GE_GTF -o $EXT_3UTR_DIR/dome.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/50pc_epiboly.stringtie.gtf     -r $GE_GTF -o $EXT_3UTR_DIR/50pc_epiboly.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/shield.stringtie.gtf           -r $GE_GTF -o $EXT_3UTR_DIR/shield.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/75pc_epiboly.stringtie.gtf     -r $GE_GTF -o $EXT_3UTR_DIR/75pc_epiboly.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/1_4_somites.stringtie.gtf      -r $GE_GTF -o $EXT_3UTR_DIR/1_4_somites.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/14_19_somites.stringtie.gtf    -r $GE_GTF -o $EXT_3UTR_DIR/14_19_somites.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/20_25_somites.stringtie.gtf    -r $GE_GTF -o $EXT_3UTR_DIR/20_25_somites.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/prim_5.stringtie.gtf           -r $GE_GTF -o $EXT_3UTR_DIR/prim_5.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/prim_15.stringtie.gtf          -r $GE_GTF -o $EXT_3UTR_DIR/prim_15.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/prim_25.stringtie.gtf          -r $GE_GTF -o $EXT_3UTR_DIR/prim_25.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/long_pec.stringtie.gtf         -r $GE_GTF -o $EXT_3UTR_DIR/long_pec.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/protruding_mouth.stringtie.gtf -r $GE_GTF -o $EXT_3UTR_DIR/protruding_mouth.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/day_4.stringtie.gtf            -r $GE_GTF -o $EXT_3UTR_DIR/day_4.ext_3UTR.bed
Rscript ./scripts/get_extended_3UTR.R -s $STR_DIR/day_5.stringtie.gtf            -r $GE_GTF -o $EXT_3UTR_DIR/day_5.ext_3UTR.bed

# Get mean and max CPM for each pervasive transcription region
EXT_3UTR_DIR="./data/stringtie/ext_3UTR"
BW_DIR="./data/BWs"
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/1_cell.ext_3UTR.bed           -f $BW_DIR/1_cell.fwd.CPM.bw           -r $BW_DIR/1_cell.rev.CPM.bw           -o $EXT_3UTR_DIR/1_cell.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/2_cell.ext_3UTR.bed           -f $BW_DIR/2_cell.fwd.CPM.bw           -r $BW_DIR/2_cell.rev.CPM.bw           -o $EXT_3UTR_DIR/2_cell.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/128_cell.ext_3UTR.bed         -f $BW_DIR/128_cell.fwd.CPM.bw         -r $BW_DIR/128_cell.rev.CPM.bw         -o $EXT_3UTR_DIR/128_cell.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/1k_cell.ext_3UTR.bed          -f $BW_DIR/1k_cell.fwd.CPM.bw          -r $BW_DIR/1k_cell.rev.CPM.bw          -o $EXT_3UTR_DIR/1k_cell.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/dome.ext_3UTR.bed             -f $BW_DIR/dome.fwd.CPM.bw             -r $BW_DIR/dome.rev.CPM.bw             -o $EXT_3UTR_DIR/dome.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/50pc_epiboly.ext_3UTR.bed     -f $BW_DIR/50pc_epiboly.fwd.CPM.bw     -r $BW_DIR/50pc_epiboly.rev.CPM.bw     -o $EXT_3UTR_DIR/50pc_epiboly.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/shield.ext_3UTR.bed           -f $BW_DIR/shield.fwd.CPM.bw           -r $BW_DIR/shield.rev.CPM.bw           -o $EXT_3UTR_DIR/shield.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/75pc_epiboly.ext_3UTR.bed     -f $BW_DIR/75pc_epiboly.fwd.CPM.bw     -r $BW_DIR/75pc_epiboly.rev.CPM.bw     -o $EXT_3UTR_DIR/75pc_epiboly.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/1_4_somites.ext_3UTR.bed      -f $BW_DIR/1_4_somites.fwd.CPM.bw      -r $BW_DIR/1_4_somites.rev.CPM.bw      -o $EXT_3UTR_DIR/1_4_somites.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/14_19_somites.ext_3UTR.bed    -f $BW_DIR/14_19_somites.fwd.CPM.bw    -r $BW_DIR/14_19_somites.rev.CPM.bw    -o $EXT_3UTR_DIR/14_19_somites.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/20_25_somites.ext_3UTR.bed    -f $BW_DIR/20_25_somites.fwd.CPM.bw    -r $BW_DIR/20_25_somites.rev.CPM.bw    -o $EXT_3UTR_DIR/20_25_somites.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/prim_5.ext_3UTR.bed           -f $BW_DIR/prim_5.fwd.CPM.bw           -r $BW_DIR/prim_5.rev.CPM.bw           -o $EXT_3UTR_DIR/prim_5.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/prim_15.ext_3UTR.bed          -f $BW_DIR/prim_15.fwd.CPM.bw          -r $BW_DIR/prim_15.rev.CPM.bw          -o $EXT_3UTR_DIR/prim_15.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/prim_25.ext_3UTR.bed          -f $BW_DIR/prim_25.fwd.CPM.bw          -r $BW_DIR/prim_25.rev.CPM.bw          -o $EXT_3UTR_DIR/prim_25.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/long_pec.ext_3UTR.bed         -f $BW_DIR/long_pec.fwd.CPM.bw         -r $BW_DIR/long_pec.rev.CPM.bw         -o $EXT_3UTR_DIR/long_pec.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/protruding_mouth.ext_3UTR.bed -f $BW_DIR/protruding_mouth.fwd.CPM.bw -r $BW_DIR/protruding_mouth.rev.CPM.bw -o $EXT_3UTR_DIR/protruding_mouth.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/day_4.ext_3UTR.bed            -f $BW_DIR/day_4.fwd.CPM.bw            -r $BW_DIR/day_4.rev.CPM.bw            -o $EXT_3UTR_DIR/day_4.ext_3UTR.CPM.bed
python ./scripts/get_extended_3UTR_CPM.py -i $EXT_3UTR_DIR/day_5.ext_3UTR.bed            -f $BW_DIR/day_5.fwd.CPM.bw            -r $BW_DIR/day_5.rev.CPM.bw            -o $EXT_3UTR_DIR/day_5.ext_3UTR.CPM.bed


# Combine extended 3' UTRs into a consensus set
#20_05_14_White2017_Find_pervasive_transcription/TES_perv_trans_regions_plots.R
EXT_3UTR_DIR="./data/stringtie/ext_3UTR"
TE_GTF_DFRAG="./data/danRer11.TEtrans_uID.dfrag.gtf"
Rscript ./scripts/combine_ext_3UTRs.R --te_gtf $TE_GTF_DFRAG \
                                      --dir_ext_UTR $EXT_3UTR_DIR \
                                      --output ${TE_GTF_DFRAG/.gtf/.ext_tmp.gtf} 



### Make TE categories
# Subset protein coding genes from gene annotations
GE_GTF="./data/Danio_rerio.GRCz11.98.chr.ERCC.addChr.gtf.gz"
cat <(grep "#\!" $GE_GTF) \
    <(grep 'gene_biotype "protein_coding"' $GE_GTF) \
    > ${GE_GTF/.gtf/.PCg.gtf}
# Run ChIPseeker to annotate TE position with respect to genes
Rscript ./scripts/run_ChIPseeker_TE_annot.R --te_gtf $TE_GTF_DFRAG \
                                            --gene_gtf ${GE_GTF/.gtf/.PCg.gtf} \
                                            --output ${TE_GTF_DFRAG/.gtf/.annot.gtf}
# Make final categories
TE_GTF_DFRAG="./data/danRer11.TEtrans_uID.dfrag.gtf"
Rscript ./scripts/make_TE_categories.R --te_gtf ${TE_GTF_DFRAG/.gtf/.annot.gtf} \
                                       --te_gtf_ext ${TE_GTF_DFRAG/.gtf/.ext_tmp.gtf} \
                                       --gene_gtf ${GE_GTF/.gtf/.PCg.gtf} \
                                       --output ${TE_GTF_DFRAG/.gtf/.classified.gtf}




### Merge telescope counts
# Join the TE counts (telescope) based on the 
# defrag_transcript_id of the defrag TE annotations.
TE_GTF_DFRAG="./data/danRer11.TEtrans_uID.dfrag.classified.gtf"
#
# fwd - fwd
python ./scripts/merge_Teles_counts_defrag.py --te_gtf $TE_GTF_DFRAG \
                                              --teles_join ./data/telescope/fwd_fwd/Edited_telescope_report/join_telesRep.final_count.tab \
                                              --output ./data/telescope/fwd_fwd/Edited_telescope_report/join_telesRep.final_count.dfrag.tab \
                                              --verbose False
# rev - rev
python ./scripts/merge_Teles_counts_defrag.py --te_gtf $TE_GTF_DFRAG \
                                              --teles_join ./data/telescope/rev_rev/Edited_telescope_report/join_telesRep.final_count.tab \
                                              --output ./data/telescope/rev_rev/Edited_telescope_report/join_telesRep.final_count.dfrag.tab \
                                              --verbose False

##################