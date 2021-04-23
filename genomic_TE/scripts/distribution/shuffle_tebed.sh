#!/usr/bin/env bash

tebed="../data/repeatmasker-out/defragmented/danRer11.nonalt.fa.out.defrag.bed"
genebed="../data/genome/Danio_rerio.GRCz11.genes.bed"
chroms="../data/genome/Danio_rerio.GRCz11.karyotype.txt"
gaps="../data/genome/Danio_rerio.GRCz11.gapfile.txt"

datadir="../data/shuffled-beds"

npermutations=50
maxprocs=6

shuffledist () {
    bedtools shuffle -i $tebed -g $chroms -excl $gaps > tmp_${1}.bed
    bedtools sort -i tmp_${1}.bed > tmp_${1}_1.bed
    mv tmp_${1}_1.bed tmp_${1}.bed
    bedtools closest -a tmp_${1}.bed -b $genebed -d \
        | awk '{ OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $NF }' \
        > ${datadir}/te_gene_shuffled_${1}.bed
    rm tmp_${i}.bed
}

# Calculate true distances
bedtools closest -a $tebed -b $genebed -d \
    | awk '{ OFS="\t" } { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $NF }' \
    > ${datadir}/te_gene_observed.bed

# Calculate shuffled distances
for i in $(seq 1 $npermutations); do
    # Limit number of processes running in background
    if [ $(jobs -r | wc -l) -ge $maxprocs ]; then
        wait $(jobs -r -p | head -n 1)
    fi

    shuffledist $i &
done

echo "Succesfully ran ${npermutations}"
