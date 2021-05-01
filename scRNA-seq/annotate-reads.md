Code for creating an expression matrix
================

**\#this is the code for the analysis of singel cell sequencing in the
manuscript entitled “A genomic portrait of zebrafish transposable
elements and their spatiotemporal embryonic expression”** <br>
\#\#re-align reads to dr11 (all file downloaded from Farrell et al., are
.bam files)

``` bash
#!/usr/bin/env bash
#/workdir/nc499/snRNA/schier/dr11/

#export JAVA_HOME=/usr/local/jdk1.8.0_121
#export PATH=$JAVA_HOME/bin:$PATH

for prefix in `cat sample.txt`; do

    # Specify filenames
    bamfile="${prefix}.bam"
    fqfile="${prefix}.fastq"
    realigned="${prefix}.realigned.sam"
    sorted=$"${prefix}.realigned.sorted.bam"
    unaligned="${prefix}.unaligned.bam"
    barcoded="${prefix}.realigned.barcoded.bam"

    echo realigning $prefix

    # Extract reads from aligned bam file.
    java -Xmx4g -jar /workdir/nc499/raw/picard/picard.jar SamToFastq \
        INPUT=$bamfile \
        FASTQ=$fqfile

    echo mapping $prefix

    # Map reads using bowtie2. First build an index, then map the reads.
    bowtie2 -x /workdir/nc499/reference/danRer11.nonalt.bowtie2idx -U $fqfile --phred33 --reorder -S $realigned

    echo Create unaligned $prefix
    # Create unaligned bam file from original. I'm not 100% sure what all of these
    # arguments do, but the important ones are REMOVE_ALIGNMENT_INFORMATION and
    # ATTRIBUTE_TO_CLEAR, especially for the XC and XM attributes.
    java -Xmx4g -jar /workdir/nc499/raw/picard/picard.jar RevertSam \
        I=$bamfile \
        O=$unaligned \
        SANITIZE=true \
        MAX_DISCARD_FRACTION=0.005 \
        ATTRIBUTE_TO_CLEAR=XT \
        ATTRIBUTE_TO_CLEAR=GE \
        ATTRIBUTE_TO_CLEAR=GS \
        ATTRIBUTE_TO_CLEAR=XF \
        ATTRIBUTE_TO_CLEAR=ZP \
        ATTRIBUTE_TO_CLEAR=XN \
        ATTRIBUTE_TO_CLEAR=AS \
        ATTRIBUTE_TO_CLEAR=OC \
        ATTRIBUTE_TO_CLEAR=OP \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=true \
        REMOVE_DUPLICATE_INFORMATION=true \
        REMOVE_ALIGNMENT_INFORMATION=true

    #create a dictionary
    #java -jar /workdir/nc499/snRNA/schier/picard/picard.jar CreateSequenceDictionary \
    #   R= danRer11.nonalt.fasta \
    #   O= danRer11.nonalt.dict

    echo sorting $prefix
    #Sort sam to bam
    java -Xmx4g -jar /workdir/nc499/raw/picard/picard.jar SortSam \
        I= $realigned \
        O= $sorted \
        SO=queryname

    echo merging $prefix
    # Merge the new bowetie2-alignment bam file with the unaligned bam file.
    java -Xmx4g -jar /workdir/nc499/raw/picard/picard.jar MergeBamAlignment \
        REFERENCE_SEQUENCE= danRer11.nonalt.fasta \
        UNMAPPED_BAM=$unaligned \
        ALIGNED_BAM=$sorted \
        OUTPUT=$barcoded \
        INCLUDE_SECONDARY_ALIGNMENTS=false \
        PAIRED_RUN=false \
        VALIDATION_STRINGENCY=SILENT
done
```

\#create a reference file for TE and gene annotation \#grep
non-readthrough loci from Quirzes gtf file

``` bash
grep 'Non-Readtrhough' TEloci_table_Cluster_tau_readtrhough.gtf > TEloci_table_Cluster_tau_NONreadtrhough.gtf
```

\#some NRT loci are double counted; clean this up

``` bash
#one NRT loci matches 2 in ref file
bedtools intersect -a <(sort -k1,1 -k2,2n TEloci_table_Cluster_tau_NONreadtrhough.gtf) -b <(sort -k1,1 -k2,2n danRer11.nonalt.fa.out.reformatted.bed) -f 0.9.8 -r -wa -c | awk '{if ($27 > 1) print}' | cut -f9 | awk '{print $4}' | sed 's/"/\\"/g' > doublecounting_RT_0729.txt

#get rid of the duplicate NRT loci
grep -v -f doublecounting_RT_0729.txt TEloci_table_Cluster_tau_NONreadtrhough.gtf > TEloci_table_Cluster_tau_NONreadtrhough_nodbcounting.gtf

#grep NRT loci from danRer11.nonalt.fa.out.reformatted.bed 
bedtools intersect -a <(sort -k1,1 -k2,2n danRer11.nonalt.fa.out.reformatted.bed) -b <(sort -k1,1 -k2,2n TEloci_table_Cluster_tau_NONreadtrhough_nodbcounting.gtf) -f 0.9.8 -r -wa > TEloci_NONrt_0729.bed
```

\#clean the format

``` bash
 #add gene_name and transcript_name---this works
awk '{print $10,$11,$12,$13}' TEloci_NONrt_0729.bed | sed 's/gene_id/gene_name/' | sed
's/transcript_id/transcript_name/' > TE_name_temp
paste -d'\ ' TEloci_NONrt_0729.bed TE_name_temp > danRer11.TEloci_NONrt_0729.genename.bed
rm TE_name_temp
#convert to gtf
awk '{print $1"\t"$7"\t"$8"\t"($2+1)"\t"$3"\t"$5"\t"$6"\t"$9"\t"(substr($0, index($0,$10)))}'
danRer11.TEloci_NONrt_0729.genename.bed  > danRer11.TEloci_NONrt_0729.genename.gtf
--------------
#change TE to exon ($3)
sed -i -r 's/\tsimilarity/\texon/' danRer11.TEloci_NONrt_0729.genename.gtf
#change TE name to individual TE
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9" "$12" "$11" "$12" "$13" "$14" "$15" "$16" "$17" "$12" "$19" "$20}'
OFS='\t'  danRer11.TEloci_NONrt_0729.genename.gtf > danRer11.TEloci_NONrt_0729.genename.dup.gtf
#combine TE and gene
cat Danio_rerio.GRCz11.intactgenename.gtf danRer11.TEloci_NONrt_0729.genename.dup.gtf | sort -k1,1 -k2,2n >
danRer11.gene.TEloci_NONrt_0729.gtf &
#this is the final file
sed 's/chrMT/chrM/g' danRer11.gene.TEloci_NONrt_0729.gtf > danRer11.gene.TEloci_NONrt.chrM.0729.gtf
sort -k1,1 -k4,4n danRer11.gene.TEloci_NONrt.chrM.0729.gtf > sorted_danRer11.gene.TEloci_NONrt.chrM.0729.gtf
rm danRer11.TEloci_NONrt_0729.genename.dup.gtf danRer11.gene.TEloci_NONrt_0729.gtf
```

\#Use dro-seq tool to create a reference file (reFLAT format) for
annotation

``` bash
/workdir/nc499/tools/Drop-seq_tools-2.3.0/ConvertToRefFlat -m 64g ANNOTATIONS_FILE= sorted_danRer11.gene.TEloci_NONrt.chrM.0729.gtf SEQUENCE_DICTIONARY= ../danRer11.nonalt.dic OUTPUT= danRer11.gene.TEloci_NONrt.chrM.0729.refFlat
```

\#re-annotate all the reads

``` bash
/workdir/nc499/tools/Drop-seq_tools-2.3.0/TagReadWithGeneFunction -m 120g INPUT=ZF75-DS2.realigned.barcoded.bam OUTPUT=ZF75-DS2.realigned.geneTE.0729.bam ANNOTATIONS_FILE=/workdir/nc499/reference/dr11_ref/danRer11.gene.TEloci_NONrt.chrM.0729.refFlat  &

/workdir/nc499/tools/Drop-seq_tools-2.3.0/TagReadWithGeneFunction -m 120g INPUT=ZF75-DS3.realigned.barcoded.bam OUTPUT=ZF75-DS3.realigned.geneTE.0729.bam ANNOTATIONS_FILE=/workdir/nc499/reference/dr11_ref/danRer11.gene.TEloci_NONrt.chrM.0729.refFlat  &

/workdir/nc499/tools/Drop-seq_tools-2.3.0/TagReadWithGeneFunction -m 120g INPUT=ZF75-DS4.realigned.barcoded.bam OUTPUT=ZF75-DS4.realigned.geneTE.0729.bam ANNOTATIONS_FILE=/workdir/nc499/reference/dr11_ref/danRer11.gene.TEloci_NONrt.chrM.0729.refFlat  &

#for loop for this
for sample in `cat sample.txt`; do /workdir/nc499/tools/Drop-seq_tools-2.3.0/TagReadWithGeneFunction -m 120g INPUT=$sample.realigned.barcoded.bam OUTPUT=$sample.realigned.geneTE.0729.bam ANNOTATIONS_FILE=/workdir/nc499/reference/dr11_ref/danRer11.gene.TEloci_NONrt.chrM.0729.refFlat; done & 
```

\#remove name (dup) for TEs, in order to combine them into one single
family

``` bash
#for loop
#replace.sam.sh
#! /bin/bash
for sample in `cat sample.txt`; 
do samtools view -h -o $sample.sam $sample;
done

#edit genebarcode.bam.firsthalf.txt to add sam
for sample in `cat sample.txt`; 
do sed -i 's/_dup[0-9]*//g' $sample.sam;
done 

#compress back to bam
for sample in `cat sample.txt`; 
do samtools view -Sb $sample.sam > ${sample/.sam/.bam};
done 
```

\#Digital expression

``` bash
#! /bin/bash
#ZF30
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF30-DS5.realigned.geneTE.0729.bam O= ZF30-DS5.dge.NRT.txt.gz SUMMARY=ZF30-DS5.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=625 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF30-DS5b.realigned.geneTE.0729.bam O= ZF30-DS5b.dge.NRT.txt.gz SUMMARY=ZF30-DS5b.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=625 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF3S
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF3S-DS5.realigned.geneTE.0729.bam O= ZF3S-DS5.dge.NRT.txt.gz SUMMARY=ZF3S-DS5.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF50
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF50-DS2.realigned.geneTE.0729.bam O= ZF50-DS2.dge.NRT.txt.gz \
SUMMARY=ZF50-DS2.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF50-DS3.realigned.geneTE.0729.bam O= ZF50-DS3.dge.NRT.txt.gz \
SUMMARY=ZF50-DS3.summary.NRT.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF50-DS4.realigned.geneTE.0729.bam O= ZF50-DS4.dge.NRT.txt.gz \
SUMMARY=ZF50-DS4.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF50-DS4b.realigned.geneTE.0729.bam O= ZF50-DS4b.dge.NRT.txt.gz \
SUMMARY=ZF50-DS4b.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF60
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF60-DS2.realigned.geneTE.0729.bam O= ZF60-DS2.dge.NRT.txt.gz \
SUMMARY=ZF60-DS2.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF60-DS3.realigned.geneTE.0729.bam O= ZF60-DS3.dge.NRT.txt.gz \
SUMMARY=ZF60-DS3.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF60-DS4.realigned.geneTE.0729.bam O= ZF60-DS4.dge.NRT.txt.gz \
SUMMARY=ZF60-DS4.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF6S
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF6S-DS5.realigned.geneTE.0729.bam O= ZF6S-DS5.dge.NRT.txt.gz \
SUMMARY=ZF6S-DS5.summary.NRT.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF6S-DS5b.realigned.geneTE.0729.bam O= ZF6S-DS5b.dge.NRT.txt.gz \
SUMMARY=ZF6S-DS5b.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF75
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF75-DS2.realigned.geneTE.0729.bam O= ZF75-DS2.dge.NRT.txt.gz \
SUMMARY=ZF75-DS2.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1400 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF75-DS3.realigned.geneTE.0729.bam O= ZF75-DS3.dge.NRT.txt.gz \
SUMMARY=ZF75-DS3.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1400 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF75-DS4.realigned.geneTE.0729.bam O= ZF75-DS4.dge.NRT.txt.gz \
SUMMARY=ZF75-DS4.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1400 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF90
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF90-DS2.realigned.geneTE.0729.bam O= ZF90-DS2.dge.NRT.txt.gz \
SUMMARY=ZF90-DS2.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF90-DS3.realigned.geneTE.0729.bam O= ZF90-DS3.dge.NRT.txt.gz \
SUMMARY=ZF90-DS3.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF90-DS4.realigned.geneTE.0729.bam O= ZF90-DS4.dge.NRT.txt.gz \
SUMMARY=ZF90-DS4.NRT.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFB
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFB-DS2.realigned.geneTE.0729.bam O= ZFB-DS2.dge.NRT.txt.gz SUMMARY=ZFB-DS2.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE 

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFB-DS2b.realigned.geneTE.0729.bam O= ZFB-DS2b.dge.NRT.txt.gz SUMMARY=ZFB-DS2b.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFB-DS3.realigned.geneTE.0729.bam O= ZFB-DS3.dge.NRT.txt.gz SUMMARY=ZFB-DS3.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFB-DS4.realigned.geneTE.0729.bam O= ZFB-DS4.dge.NRT.txt.gz SUMMARY=ZFB-DS4.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFDOME
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFDOME-DS5.realigned.geneTE.0729.bam O= ZFDOME-DS5.dge.NRT.txt.gz SUMMARY=ZFDOME-DS5.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=800 MIN_NUM_READS_PER_CELL=2000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFHIGH
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFHIGH-DS5.realigned.geneTE.0729.bam O= ZFHIGH-DS5.dge.NRT.txt.gz SUMMARY=ZFHIGH-DS5.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=1000 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFHIGH-DS5b.realigned.geneTE.0729.bam O= ZFHIGH-DS5b.dge.NRT.txt.gz SUMMARY=ZFHIGH-DS5b.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=1000 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFOBLONG
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFOBLONG-DS5.realigned.geneTE.0729.bam O= ZFOBLONG-DS5.dge.NRT.txt.gz SUMMARY=ZFOBLONG-DS5.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=625 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFOBLONG-DS5b.realigned.geneTE.0729.bam O= ZFOBLONG-DS5b.dge.NRT.txt.gz SUMMARY=ZFOBLONG-DS5b.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=625 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFS
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFS-DS5.realigned.geneTE.0729.bam O= ZFS-DS5.dge.NRT.txt.gz SUMMARY=ZFS-DS5.NRT.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE
```

\#replace the title

``` bash
#! /bin/bash
sed -i -e '1s/\t/\tZF30-DS5_/g' ZF30-DS5.dge.NRT.txt 
sed -i -e '1s/\t/\tZF30-DS5b_/g'    ZF30-DS5b.dge.NRT.txt
sed -i -e '1s/\t/\tZF3S-DS5_/g' ZF3S-DS5.dge.NRT.txt 
sed -i -e '1s/\t/\tZF50-DS2_/g' ZF50-DS2.dge.NRT.txt 
sed -i -e '1s/\t/\tZF50-DS3_/g' ZF50-DS3.dge.NRT.txt 
sed -i -e '1s/\t/\tZF50-DS4_/g' ZF50-DS4.dge.NRT.txt 
sed -i -e '1s/\t/\tZF50-DS4b_/g'    ZF50-DS4b.dge.NRT.txt
sed -i -e '1s/\t/\tZF60-DS2_/g' ZF60-DS2.dge.NRT.txt 
sed -i -e '1s/\t/\tZF60-DS3_/g' ZF60-DS3.dge.NRT.txt 
sed -i -e '1s/\t/\tZF60-DS4_/g' ZF60-DS4.dge.NRT.txt 
sed -i -e '1s/\t/\tZF6S-DS5_/g' ZF6S-DS5.dge.NRT.txt 
sed -i -e '1s/\t/\tZF6S-DS5b_/g'    ZF6S-DS5b.dge.NRT.txt
sed -i -e '1s/\t/\tZF75-DS2_/g' ZF75-DS2.dge.NRT.txt 
sed -i -e '1s/\t/\tZF75-DS3_/g' ZF75-DS3.dge.NRT.txt 
sed -i -e '1s/\t/\tZF75-DS4_/g' ZF75-DS4.dge.NRT.txt 
sed -i -e '1s/\t/\tZF90-DS2_/g' ZF90-DS2.dge.NRT.txt 
sed -i -e '1s/\t/\tZF90-DS3_/g' ZF90-DS3.dge.NRT.txt 
sed -i -e '1s/\t/\tZF90-DS4_/g' ZF90-DS4.dge.NRT.txt 
sed -i -e '1s/\t/\tZFB-DS2_/g'  ZFB-DS2.dge.NRT.txt
sed -i -e '1s/\t/\tZFB-DS2b_/g' ZFB-DS2b.dge.NRT.txt
sed -i -e '1s/\t/\tZFB-DS3_/g'  ZFB-DS3.dge.NRT.txt
sed -i -e '1s/\t/\tZFB-DS4_/g'  ZFB-DS4.dge.NRT.txt
sed -i -e '1s/\t/\tZFDOME-DS5_/g'   ZFDOME-DS5.dge.NRT.txt
sed -i -e '1s/\t/\tZFHIGH-DS5_/g'   ZFHIGH-DS5.dge.NRT.txt
sed -i -e '1s/\t/\tZFHIGH-DS5b_/g'  ZFHIGH-DS5b.dge.NRT.txt
sed -i -e '1s/\t/\tZFOBLONG-DS5_/g' ZFOBLONG-DS5.dge.NRT.txt
sed -i -e '1s/\t/\tZFOBLONG-DS5b_/g'    ZFOBLONG-DS5b.dge.NRT.txt
sed -i -e '1s/\t/\tZFS-DS5_/g'  ZFS-DS5.dge.NRT.txt
```

\#filter cells in R, to get rid of cells with abnormal (too low/high)
transcipt numbers

``` bash
library(Seurat)
library(dplyr)
library(data.table)

data.orig <- read.csv("sub_ZFS.txt", sep = "\t", row.names = 1, header = T )  #substitue the sample name for each sample
data.obj <- CreateSeuratObject(counts = data.orig, project ="filter")
data.obj
data.obj<-subset(data.obj, subset = nFeature_RNA > 600 & nFeature_RNA < 2700)
data.obj
data_to_write_out <- as.data.frame(as.matrix(data.obj[["RNA"]]@data))
fwrite(x = data_to_write_out, sep="\t", row.names = TRUE, file = "sub_ZFS.clean.matrix")
rm(data.orig,data.obj,data_to_write_out)
```

\#Jon helped me to merge all samples into one big matrix \#named
filtered\_clean.txt

``` bash
sed -i 's/nan/0.0/g' filtered_clean.txt
grep 'ENSDAR' filtered_clean.txt | awk '{ for(i=2; i<=NF;i++) j+=$i; print j; j=0 }' | awk '{total += $1; count++ } END {print total}'
255012447
grep -v 'ENSDAR' filtered_clean.txt | awk '{ for(i=2; i<=NF;i++) j+=$i; print j; j=0 }' | awk '{total += $1; count++ } END {print total}'
5613126
2.2%
```
