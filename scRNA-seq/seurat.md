Code for identifying cell-specific TEs (Seurat)
================

**\#this is the code for the analysis of singel cell sequencing in the
manuscript entitled “A genomic portrait of zebrafish transposable
elements and their spatiotemporal embryonic expression”** <br>
\#\#create cell meta

``` bash
head -1 filtered_clean_new_annot2.txt | sed 's/\t/\n/g' | awk '{print $1,$1}' | sed 's/\./\t/2' | sed 's/_/\t/2' | awk -v OFS='\t' '{print $1,$3}' > cell_meta.txt
```

\#\#Seurat

``` r
library(dplyr)
library(Seurat)
setwd("/workdir/nc499/raw/new_annot/dge/filtered")
pbmc.data <- read.csv("filtered_clean_new_annot2.txt", sep = "\t", row.names = 1, header = T )
pbmc <- CreateSeuratObject(pbmc.data, project = "all", min.cells = 3, min.features = 100)#change the min.features
dim(pbmc) 


#remove MT
MT<-read.delim("MT_genes_0130.txt", header=TRUE)
all(MT %in% rownames(pbmc))

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = MT$genenames)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 100 & percent.mt < 45 )
dim(pbmc)

#add stages
stages=read.table("cell_stages.txt", sep="\t", header=FALSE,row.names=1)
pbmc<-AddMetaData(pbmc, stages, col.name="stages")
head(pbmc@meta.data)

#add batch data
cellbatch=read.table("cell_meta.txt", sep="\t", header=FALSE,row.names=1)
pbmc<-AddMetaData(pbmc, cellbatch, col.name="batch")
head(pbmc@meta.data)

#separate batches
pbmc.list <- SplitObject(pbmc, split.by = "batch")

for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE)
  pbmc.list[[i]] <- FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", 
                                         nfeatures = 5000, verbose = TRUE)
}

#options(future.globals.maxSize = 10000 * 1024^2)

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features=5000, dims = 1:40)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:40)
DefaultAssay(pbmc.integrated) <- "integrated"
pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)


pbmc.integrated <- RunPCA(pbmc.integrated, npcs = 300)
pbmc.integrated <- JackStraw(pbmc.integrated, reduction = "pca",num.replicate =100,dims = 150) #num.replicate = 600, dims = 300

pbmc.integrated <- ScoreJackStraw(pbmc.integrated, dims = 1:150)

JackStrawPlot(pbmc.integrated, dims = 1:150)

ElbowPlot(pbmc.integrated, ndims = 150)
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:143) 
pbmc.integrated <- FindClusters(pbmc.integrated, resolution =5.8) #default= resolution =0.5, old 5.6
#head(Idents(pbmc), 5)

pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:143)

DimPlot(pbmc.integrated, reduction = "umap", label =T,label.size=3,group.by="orig.ident", cols = viridis(28))#group.by="orig.ident",, cols = viridis(28)

FeaturePlot(pbmc.integrated, features = c("ENSDARG00000068255","ENSDARG00000022813","ENSDARG00000014373"))
FeaturePlot(pbmc.integrated, features = c(
  "BHIKHARI-3-LTR-DR","BHIKHARI-5-LTR-DR","BHIKHARI-LTR","ERV1-3-I-DR"),
  cols = c("honeydew3","red1"),pt.size=0.3,min.cutoff=0)
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

saveRDS(pbmc.integrated,file="all_dr11_dm143r5.8_batchMT_0724.rds")


#find markers
pbmc.integrated.markers.20 <- FindAllMarkers(pbmc.integrated, min.pct=0.2, min.diff.pct=0.2, logfc.threshold = 0.25, return.thresh=0.05, only.pos=TRUE) 
pbmc.markers.20.2<-pbmc.integrated.markers.20[which(pbmc.integrated.markers.20[,5]<0.05),]

write.table(pbmc.markers.20.2,file="all_r5.8_dim143_batchcorrectedMT_markers0724_0.2_pos_0.2diff.txt",sep="\t",quote=F,row.names=F)

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

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:50)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution =1.2) 
pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:50)
DimPlot(pbmc, reduction = "umap", label =T, label.size=4)#,group.by="orig.ident"
FeaturePlot(pbmc.integrated, features = c("ENSDARG00000068255","ENSDARG00000022813","ENSDARG00000014373"))

pbmc.list
AE<-AverageExpression(object = pbmc.integrated)
apply(AE$RNA,2,sum)
write.table(AE$RNA,file="blastula_r1.2_dm50_batchcorrected_AE.txt",sep="\t",quote=F,row.names=T)


##gastrula
nfeatures = 3000, anchor dims = 1:30
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:90)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution =3) 

##segmentation
nfeatures = 3000, anchor dims = 1:30, dim60 r2
```
