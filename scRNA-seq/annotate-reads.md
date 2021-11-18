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
grep 'Self_expressed' danRer11.TEtrans_uID.subClass.wredundant.annot_categ.PCg.dfrag.c2.Rt.cluster.gtf > danRer11.TEtrans_uID.subClass_SelfE.0706.gtf
```

\#combine with my old annotation
\#make the format is readable for redFLAT

``` {bash, eval=FALSE}
sed 's/$/;/' danRer11.TEtrans_uID.subClass_SelfE.0706.gtf > danRer11.TEtrans_uID.subClass_SelfE.0706.2.gtf
awk '{print $9,$10,$11,$12}' danRer11.TEtrans_uID.subClass_SelfE.0706.2.gtf | sed 's/gene_id/gene_name/' | sed 's/transcript_id/transcript_name/' > TE_name_temp
paste -d'\ ' danRer11.TEtrans_uID.subClass_SelfE.0706.2.gtf TE_name_temp > danRer11.TEtrans_uID.subClass_SelfE.genename.0706.gtf
rm TE_name_temp
rm danRer11.TEtrans_uID.subClass_SelfE.0706.2.gtf

#change TE name to individual TE
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9" "$12" "$11" "$12" "$13" "$14" "$15" "$16" "$17" "$18" "$55" "$12" "$57" "$58}' OFS='\t'  danRer11.TEtrans_uID.subClass_SelfE.genename.0706.gtf 
> danRer11.TEtrans_uID.subClass_SelfE.indTE.0706.gtf

#combine TE and gene 
cat /workdir/nc499/reference/dr11_ref/Danio_rerio.GRCz11.intactgenename.gtf danRer11.TEtrans_uID.subClass_SelfE.indTE.0706.gtf | sort -k1,1 -k2,2n > danRer11.gene.TESelfE_0706.gtf

#this is the final file
sed -i 's/chrMT/chrM/g' danRer11.gene.TESelfE_0706.gtf

#Use dro-seq tool to create a reference file (reFLAT format) for annotation
/workdir/nc499/tools/Drop-seq_tools-2.3.0/ConvertToRefFlat -m 64g ANNOTATIONS_FILE= danRer11.gene.TESelfE_0706.gtf SEQUENCE_DICTIONARY= ../danRer11.nonalt.dic OUTPUT= danRer11.gene.TESelfE_0706.refFlat 

```

#re-annotate all the reads
```{bash, eval=FALSE}
for sample in `cat sample.txt`; do 
/workdir/nc499/tools/Drop-seq_tools-2.3.0/TagReadWithGeneFunction \
-m 120g \
INPUT=../$sample \
OUTPUT=$sample.realigned.geneTEbarcoded.070821.bam \
ANNOTATIONS_FILE=/workdir/nc499/reference/dr11_ref/danRer11.gene.TESelfE_0706.refFlat; 
done & 

```

#remove name (dup) for TEs, in order to combine them into one single family
```{bash, eval=FALSE}
#for loop

#replace.sam.sh
#! /bin/bash/env
for sample in `cat sample.txt`; 
do samtools view -h $sample.realigned.geneTEbarcoded.070821.bam | sed 's/_dup[0-9]*//g' | samtools view -Sb - > $sample.family.bam;
done
```

\#Digital expression

``` {bash, eval=FALSE}
#! /bin/bash
#ZF30  
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF30-DS5.realigned.barcoded.bam.family.bam O= ZF30-DS5.dge.newannot.txt.gz \
SUMMARY=ZF30-DS5.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=625 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF30-DS5b.realigned.barcoded.bam.family.bam O= ZF30-DS5b.dge.newannot.txt.gz \
SUMMARY=ZF30-DS5b.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=625 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF3S
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF3S-DS5.realigned.barcoded.bam.family.bam O= ZF3S-DS5.dge.newannot.txt.gz \
SUMMARY=ZF3S-DS5.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF50
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF50-DS2.realigned.barcoded.bam.family.bam O= ZF50-DS2.dge.newannot.txt.gz \
SUMMARY=ZF50-DS2.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF50-DS3.realigned.barcoded.bam.family.bam O= ZF50-DS3.dge.newannot.txt.gz \
SUMMARY=ZF50-DS3.summary.newannot.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF50-DS4.realigned.barcoded.bam.family.bam O= ZF50-DS4.dge.newannot.txt.gz \
SUMMARY=ZF50-DS4.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF50-DS4b.realigned.barcoded.bam.family.bam O= ZF50-DS4b.dge.newannot.txt.gz \
SUMMARY=ZF50-DS4b.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF60
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF60-DS2.realigned.barcoded.bam.family.bam O= ZF60-DS2.dge.newannot.txt.gz \
SUMMARY=ZF60-DS2.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF60-DS3.realigned.barcoded.bam.family.bam O= ZF60-DS3.dge.newannot.txt.gz \
SUMMARY=ZF60-DS3.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF60-DS4.realigned.barcoded.bam.family.bam O= ZF60-DS4.dge.newannot.txt.gz \
SUMMARY=ZF60-DS4.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF6S
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF6S-DS5.realigned.barcoded.bam.family.bam O= ZF6S-DS5.dge.newannot.txt.gz \
SUMMARY=ZF6S-DS5.summary.newannot.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF6S-DS5b.realigned.barcoded.bam.family.bam O= ZF6S-DS5b.dge.newannot.txt.gz \
SUMMARY=ZF6S-DS5b.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF75
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF75-DS2.realigned.barcoded.bam.family.bam O= ZF75-DS2.dge.newannot.txt.gz \
SUMMARY=ZF75-DS2.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1400 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF75-DS3.realigned.barcoded.bam.family.bam O= ZF75-DS3.dge.newannot.txt.gz \
SUMMARY=ZF75-DS3.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1400 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF75-DS4.realigned.barcoded.bam.family.bam O= ZF75-DS4.dge.newannot.txt.gz \
SUMMARY=ZF75-DS4.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1400 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZF90
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF90-DS2.realigned.barcoded.bam.family.bam O= ZF90-DS2.dge.newannot.txt.gz \
SUMMARY=ZF90-DS2.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF90-DS3.realigned.barcoded.bam.family.bam O= ZF90-DS3.dge.newannot.txt.gz \
SUMMARY=ZF90-DS3.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZF90-DS4.realigned.barcoded.bam.family.bam O= ZF90-DS4.dge.newannot.txt.gz \
SUMMARY=ZF90-DS4.newannot.summary.txt TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 \
MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#! /bin/bash

#ZFB
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFB-DS2.realigned.barcoded.bam.family.bam O= ZFB-DS2.dge.newannot.txt.gz SUMMARY=ZFB-DS2.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE 

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFB-DS2b.realigned.barcoded.bam.family.bam O= ZFB-DS2b.dge.newannot.txt.gz SUMMARY=ZFB-DS2b.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFB-DS3.realigned.barcoded.bam.family.bam O= ZFB-DS3.dge.newannot.txt.gz SUMMARY=ZFB-DS3.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFB-DS4.realigned.barcoded.bam.family.bam O= ZFB-DS4.dge.newannot.txt.gz SUMMARY=ZFB-DS4.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=500 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFDOME
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFDOME-DS5.realigned.barcoded.bam.family.bam O= ZFDOME-DS5.dge.newannot.txt.gz SUMMARY=ZFDOME-DS5.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=800 MIN_NUM_READS_PER_CELL=2000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFHIGH
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFHIGH-DS5.realigned.barcoded.bam.family.bam O= ZFHIGH-DS5.dge.newannot.txt.gz SUMMARY=ZFHIGH-DS5.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=1000 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFHIGH-DS5b.realigned.barcoded.bam.family.bam O= ZFHIGH-DS5b.dge.newannot.txt.gz SUMMARY=ZFHIGH-DS5b.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=1000 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFOBLONG
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFOBLONG-DS5.realigned.barcoded.bam.family.bam O= ZFOBLONG-DS5.dge.newannot.txt.gz SUMMARY=ZFOBLONG-DS5.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=625 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFOBLONG-DS5b.realigned.barcoded.bam.family.bam O= ZFOBLONG-DS5b.dge.newannot.txt.gz SUMMARY=ZFOBLONG-DS5b.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=625 MIN_NUM_READS_PER_CELL=1500 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE

#ZFS
/workdir/nc499/tools/Drop-seq_tools-2.3.0/DigitalExpression -m 120g \
I=ZFS-DS5.realigned.barcoded.bam.family.bam O= ZFS-DS5.dge.newannot.txt.gz SUMMARY=ZFS-DS5.newannot.summary.txt \
TMP_DIR=/workdir/nc499/temp/ READ_MQ= -1 MIN_NUM_GENES_PER_CELL=600 MIN_NUM_READS_PER_CELL=1000 USE_STRAND_INFO=false STRAND_STRATEGY=SENSE
```



\#filter cells in R, to get rid of cells with abnormal (too low/high)
transcipt numbers

``` bash
library(Seurat)
library(dplyr)
library(data.table)

data.orig <- read.csv("ZFS-DS5.dge.newannot.txt", sep = "\t", row.names = 1, header = T )  #substitue the sample name for each sample
data.obj <- CreateSeuratObject(counts = data.orig, project ="filter")
data.obj
data.obj<-subset(data.obj, subset = nFeature_RNA > 600 & nFeature_RNA < 2700)
data.obj
data_to_write_out <- as.data.frame(as.matrix(data.obj[["RNA"]]@data))
fwrite(x = data_to_write_out, sep="\t", row.names = TRUE, file = "sub_ZFS.clean.matrix")
rm(data.orig,data.obj,data_to_write_out)
```

\#Jon helped me to merge all samples into one big matrix \#named
filtered_clean_new_annot2.txt
