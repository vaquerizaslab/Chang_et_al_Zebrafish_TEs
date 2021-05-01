Code for identifying cell-specific TEs (Seurat)
================

**\#this is the code for the analysis of singel cell sequencing in the
manuscript entitled “A genomic portrait of zebrafish transposable
elements and their spatiotemporal embryonic expression”** <br>
\#\#create cell meta

``` bash
head -1 filtered_clean.txt | sed 's/\t/\n/g' | awk '{print $1,$1}' | sed 's/\./\t/2' | sed 's/_/\t/2' | awk -v OFS='\t' '{print $1,$3}' > cell_meta.txt
```

\#\#Seurat

``` r
library(dplyr)
library(Seurat)
setwd("/workdir/nc499/raw/NRT/filtered")
pbmc.data <- read.csv("filtered_clean.txt", sep = "\t", row.names = 1, header = T )
pbmc <- CreateSeuratObject(pbmc.data, project = "all", min.cells = 3, min.features = 100)#change the min.features

#remove MT
MT<-read.delim("MT_genes_0130.txt")
all(MT %in% rownames(pbmc))

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = MT$genenames)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 100 & percent.mt < 45 )

#add stages
stages=read.table("cell_stages.txt", sep="\t", header=FALSE,row.names=1)
pbmc<-AddMetaData(pbmc, stages, col.name="stages")

#add batch data
cellbatch=read.table("cell_meta.txt", sep="\t", header=TRUE,row.names=1)
pbmc<-AddMetaData(pbmc, cellbatch, col.name="batch")

#separate batches
pbmc.list <- SplitObject(pbmc, split.by = "batch")

for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE)
  pbmc.list[[i]] <- FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", 
                                         nfeatures = 5000, verbose = TRUE)
}

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features=5000, dims = 1:40)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:40)
DefaultAssay(pbmc.integrated) <- "integrated"
pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)


pbmc.integrated <- RunPCA(pbmc.integrated, npcs = 300)
pbmc.integrated <- JackStraw(pbmc.integrated, reduction = "pca",num.replicate =600,dims = 300 ) #num.replicate = 600, dims = 300

pbmc.integrated <- ScoreJackStraw(pbmc.integrated, dims = 1:150)

JackStrawPlot(pbmc.integrated, dims = 1:150)

ElbowPlot(pbmc.integrated, ndims = 150)
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:114)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution =5.6) #default= resolution =0.5
#head(Idents(pbmc), 5)

pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:114)

DimPlot(pbmc.integrated, reduction = "umap", label =T,label.size=2)#group.by="orig.ident",, cols = viridis(28)

FeaturePlot(pbmc.integrated, features = c("ENSDARG00000068255","ENSDARG00000022813","ENSDARG00000014373"))
FeaturePlot(pbmc.integrated, features = c(
  "DNA2-8-DR","DNA6-10-DR","Tx1-13-DR","Tx1-58-DR"))
#neural
FeaturePlot(pbmc.integrated, features = c(
  "ENSDARG00000003411", "ENSDARG00000038867", "ENSDARG00000068567", "ENSDARG00000101919", "ENSDARG00000043923", "ENSDARG00000077467"
))
#axial mesoderm
FeaturePlot(pbmc.integrated, features = c(
  "ENSDARG00000040944", "ENSDARG00000021201", "ENSDARG00000101576", "ENSDARG00000077403", "ENSDARG00000069093", "ENSDARG00000071082", "ENSDARG00000019122", "ENSDARG00000023656", "ENSDARG00000039173"
))

#EVL and apoptosis

"ENSDARG00000017624","ENSDARG00000058371","ENSDARG00000094041"
FeaturePlot(pbmc.integrated, features = c(
  "ENSDARG00000086374","ENSDARG00000042904","ENSDARG00000043581","ENSDARG00000017624","ENSDARG00000058371","ENSDARG00000094041"
))

#nanos3=ENSDARG00000068255, krt4=ENSDARG00000017624, noto=ENSDARG00000021201, fsta=ENSDARG00000052846
FeaturePlot(pbmc.integrated, features = c(
  "ENSDARG00000068255","ENSDARG00000017624","ENSDARG00000021201","ENSDARG00000052846"
),cols = c("honeydew3","red1"),pt.size=0.3,min.cutoff=0)



#stages
library(viridis)
pbmc.integrated@meta.data$stages <- factor(pbmc.integrated@meta.data$stages, levels=c("ZFHIGH","ZFOBLONG","ZFDOME","ZF30","ZF50","ZFS","ZF60","ZF75","ZF90","ZFB","ZF3S","ZF6S"))
DimPlot(pbmc.integrated,reduction = "umap",group.by="stages", cols = viridis(12,alpha=0.5, direction = -1))#,group.by="orig.ident",, label =T, label.size=2
#

saveRDS(pbmc.integrated,file="all_dr11_dm114r5.6_batchMT_0817.rds")


#find markers
pbmc.integrated.markers <- FindAllMarkers(pbmc.integrated, min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.05)
pbmc.markers2<-pbmc.integrated.markers[which(pbmc.integrated.markers[,5]<0.05),]

write.table(pbmc.markers2,file="all_r5.6_pc114_batchcorrectedMT_markers0817_0.1.txt",sep="\t",quote=F,row.names=F)
```

\#\#stage heatmap <br> \#\#separate stages

``` bash
perl -lane '$,="\t";
   $. == 1 and @A = grep $F[$_] =~ /ZFHIGH|ZFOBLONG|ZFDOME|ZF30/, 0..$#F;
   print @F[0,@A];
' after_filtered.matrix > blastula.matrix

perl -lane '$,="\t";
   $. == 1 and @A = grep $F[$_] =~ /ZF50|ZF60|ZFS|ZF75|ZF90|ZFB/, 0..$#F;
   print @F[0,@A];
' after_filtered.matrix > gastrula.matrix

perl -lane '$,="\t";
   $. == 1 and @A = grep $F[$_] =~ /ZF3S|ZF6S/, 0..$#F;
   print @F[0,@A];
' after_filtered.matrix > segmentation.matrix
```

\#Seurat

``` r
##blastula
pbmc.data <- read.csv("blastula.matrix", sep = "\t", row.names = 1, header = T )
pbmc <- CreateSeuratObject(pbmc.data, project = "all", min.cells = 3, min.features = 100)#change the min.features

pbmc.list <- SplitObject(pbmc, split.by = "orig.ident")
pbmc.list
for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE, scale.factor = 1000) #add scale factor
  pbmc.list[[i]] <- FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", 
                                         nfeatures = 3000, verbose = FALSE)
}


pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features=3000, dims = 1:30,k.filter=30) #k.filter=50,,anchor.features=38000,reduction="cca"
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:30)
DefaultAssay(pbmc.integrated) <- "integrated"
pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunPCA(pbmc.integrated, npcs=100)
pbmc.integrated <- JackStraw(pbmc.integrated, reduction = "pca",num.replicate =100, dim=80) 
pbmc.integrated <- ScoreJackStraw(pbmc.integrated, dims = 1:80)
JackStrawPlot(pbmc.integrated, dims = 1:80)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:43)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution =1.2) 
pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:43)
DimPlot(pbmc, reduction = "umap", label =T, label.size=4)#,group.by="orig.ident"
FeaturePlot(pbmc.integrated, features = c("ENSDARG00000068255","ENSDARG00000022813","ENSDARG00000014373"))

pbmc.list
AE<-AverageExpression(object = pbmc.integrated)
apply(AE$RNA,2,sum)
write.table(AE$RNA,file="blastula_r1.2_dm43_batchcorrected_AE.txt",sep="\t",quote=F,row.names=T)

pbmc.integrated.markers <- FindAllMarkers(pbmc.integrated, min.pct = 0.1, logfc.threshold = 0.25, return.thresh=0.05)
pbmc.markers2<-pbmc.integrated.markers[which(pbmc.integrated.markers[,5]<0.05),]

write.table(pbmc.markers2,file="blastula_r1.2_dm43_batchcorrected_markers.txt",sep="\t",quote=F,row.names=F)
saveRDS(pbmc.integrated,file="blastula_r1.2_dm43_batchcorrected_0817.rds")

##gastrula
nfeatures = 3000, anchor dims = 1:30
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:103)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution =1.5) 
42 clusters

##segmentation
nfeatures = 3000, anchor dims = 1:30, dim63 r1.5
```
