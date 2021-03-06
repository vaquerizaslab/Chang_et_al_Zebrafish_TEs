Code for the URD tree
================

**\#this is the code for the analysis of single cell sequencing in the
manuscript entitled “A genomic portrait of zebrafish transposable
elements and their spatiotemporal embryonic expression”**

\#\#Follow the supplementary analysis from Farrell et al., 2018

``` r
library(URD)
library(gridExtra)
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))
knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()
setwd("/workdir/nc499/raw/NRT/filtered")
count.all <- as.matrix(read.table("after_filtered.matrix")) 
meta.all <- read.table("cells_afterfiltered_meta.txt")
object <- createURD(count.data = count.all, meta = meta.all, min.cells=20, gene.max.cut=10000)

rm(list=c("count.all","meta.all"))
shhhh<- gc()
#find variable genes
stages<-unique(object@meta$STAGE)
cell.each.stage <- lapply(stages, function(stage) rownames(object@meta)[which(object@meta$STAGE == stage)])
var.genes.by.stage<- lapply(1:length(stages), function(n) findVariableGenes(object,cells.fit = cell.each.stage[[n]], set.object.var.genes = F, diffCV.cutoff = 0.3, mean.min= 0.005, mean.max = 100, main.use = stages[[n]], do.plot = T))
warnings()

names(var.genes.by.stage)<-stages
#take union of variable genes from all stages
var.genes<-sort(unique(unlist(var.genes.by.stage)))
#set variable genes in object
object@var.genes<-var.genes
#save variable gene lists

for (stage in stages) {
  write(var.genes.by.stage[[stage]], file=paste0("var_gene/var_", stage, ".txt"))
}
write(var.genes, file= "var_genes.txt")
#------------------------------------------------------------
object <- calcPCA(object)
set.seed(18)
object<-calcTsne(object, perplexity=30, theta=0.5)
set.seed(17)
object<-graphClustering(object, dim.use="pca", num.nn = c(15, 20, 30), do.jaccard = T, method = "Louvain")
#Need to make a new version of stage names that is alphabetical
object@meta$stage.nice <- plyr::mapvalues(x= object@meta$STAGE, from = c("ZFHIGH", "ZFOBLONG","ZFDOME","ZF30","ZF50","ZFS","ZF60", "ZF75", "ZF90", "ZFB", "ZF3S", "ZF6S"), to =c("A-HIGH", "B-OBLONG","C-DOME", "D-30", "E-50", "F-S", "G-60", "H-75", "I-90", "J-B", "K-S3", "L-6S"))

object@meta$STAGE
object@meta$stage.nice

stage.colors <- c("#CCCCCC", RColorBrewer::brewer.pal(9, "Set1")[9], RColorBrewer::brewer.pal(12, "Paired")[c(9, 10, 7, 8, 5, 6, 3, 4, 1, 2)])

plotDim(object, "stage.nice", discrete.colors=stage.colors, legend=T, plot.title= "Development Stage", alpha=0.5)
plotDim(object, "Louvain-15", legend=T, plot.title= "Louvain-Jaccard Graph-based Clustering (15 NNs)", alpha=1) 
#remove outliers
object <- calcKNN(object,nn=100)
outliers <- knnOutliers(object, nn.1=1, nn.2=20, x.max=48,slope.r=1.1, int.r=2.9, slope.b= 0.87, int.b = 10, title = "Identifying Outliers by k-NN Distance")

#rm(outliers)

gridExtra::grid.arrange(grobs=list(
  plotDim(object, "ENSDARG00000086374", alpha=0.4, point.size=0.5),
  plotDim(object, "ENSDARG00000042904", alpha=0.4, point.size=0.5),
  plotDim(object, "ENSDARG00000043581", alpha=0.4, point.size=0.5),
  plotDimHighlight(object,clustering="Louvain-15", cluster="25",legend=F)
))

#gridExtra::grid.arrange(grobs=list(
  plotDimHighlight(object,clustering="Louvain-15", cluster="27",legend=F),
  plotDimHighlight(object,clustering="Louvain-15", cluster="24",legend=F),
  plotDimHighlight(object,clustering="Louvain-15", cluster="25",legend=F),
  plotDimHighlight(object,clustering="Louvain-15", cluster="26",legend=F)
  ))




apoptotic.like.cells <- cellsInCluster(object,"Louvain-15","25")
cells.keep <- setdiff(colnames(object@logupx.data), c(outliers, apoptotic.like.cells))
object <- urdSubset(object,cells.keep=cells.keep)
saveRDS(object,file="URD_object_all_dr11_trimmed_082520.rds")
#-----------------------------------------------------------
#Diffusion map and pseudotime
object <- readRDS("URD_object_all_dr11_trimmed_082520.rds")
object<-calcDM(object,knn=200,sigma.use=11) 
saveRDS(object@dm, file="all_dr11_dm11.rds")

stage.colors <- c("#CCCCCC", RColorBrewer::brewer.pal(9, "Set1")[9], RColorBrewer::brewer.pal(12, "Paired")[c(9, 10, 7, 8, 5, 6, 3, 4, 1, 2)])

plotDimArray(object=object, reduction.use= "dm", dims.to.plot=1:11, label = "stage.nice", plot.title = "", outer.title = "STAGE -Diffusion Map Sigma 11", legend = F, alpha = 0.3, discrete.colors = stage.colors)
root.cells <-rownames(object@meta)[object@meta$STAGE == "ZFHIGH"]
#do the flood
flood.result <- floodPseudotime(object, root.cells=root.cells, n=150, minimum.cells.flooded=2, verbose=T)
#save the result
saveRDS(flood.result, file ="flood-dm-11-dr11-local.rds")

#check the result
floods.dm11 <-lapply(list.files(path="./", pattern="flood-dm-11-dr11-local",full.names=T), readRDS)
object <-floodPseudotimeProcess(object, floods.dm11, floods.name="pseudotime", max.frac.NA=0.4, pseudotime.fun=mean, stability.div=20)
pseudotimePlotStabilityOverall(object)

#inspect pseudotime
plotDim(object,"pseudotime", plot.title="Pseudotime")
#define a properyordered stage name
object@meta$HPFSTAGE <- apply(object@meta[, c("HPF", "STAGE")], 1, paste0, collapse = "-", sep= "")
#create a data.frame that includes pseudotime and stage information
gg.data <- cbind(object@pseudotime, object@meta[rownames(object@pseudotime), ])
#plot
ggplot(gg.data, aes(x = pseudotime, color = HPFSTAGE, fill = HPFSTAGE)) + geom_density(alpha = 0.4) + theme_bw()

saveRDS(object, file="object_3_withDMandPT_dm11.rds")
#####################################################################
#URD3 determining tips
library(URD)
library(scran)
library(gridExtra)
object<-readRDS("object_3_withDMandPT_dm11.rds")

#trim to cells from final stage
#subset the object
cells.6s <- grep("ZF6S", colnames(object@logupx.data), value=T)
object.6s <- urdSubset(object, cells.keep= cells.6s)
print(object.6s@logupx.data)
print(object.6s)
#perform PCA/tSNE on final stage cells
#load the variable genes specific to this stage
var.genes.6s <- scan("var_gene/var_ZF6S.txt", what= "character")

object.6s@var.genes <- var.genes.6s

print(object.6s@var.genes)
#Calculate PCA
object.6s <- calcPCA(object.6s)

#calculate tSNE
set.seed(18)
object.6s <-calcTsne(object.6s, dim.use= "pca", perplexity=30, theta=0.5)

#look at batch information
plotDim(object.6s, "BATCH", plot.title = "BATCH (Uncorrected)")

#Batch correct data
#make a copy of the object
object.6s.mnn <-object.6s
#generate expression matrices from each batch
cells.1 <-rownames(object.6s.mnn@meta)[which(object.6s.mnn@meta$BATCH == "DS5")]
cells.2 <-rownames(object.6s.mnn@meta)[which(object.6s.mnn@meta$BATCH == "DS5b")]
exp.1 <-as.matrix(object.6s.mnn@logupx.data[,cells.1])
exp.2 <-as.matrix(object.6s.mnn@logupx.data[,cells.2])

#batch correct using MNN, correcting all genes, but using the variable genes to determine mutual nearest neighbors
logupx.6s.mnn <- mnnCorrect(exp.1, exp.2, subset.row=object.6s.mnn@var.genes, k=20, sigma=1, svd.dim = 0, cos.norm.in = T, cos.norm.out=F)

#combine the resultant corrected matrices and return to original order of cells
logupx.6s.mnn<-do.call("cbind",logupx.6s.mnn[[1]])
logupx.6s.mnn<-logupx.6s.mnn[, colnames(object.6s.mnn@logupx.data)]
#re-sparsify the matrix, by turning anything less than 0 or near 0 back to 0
logupx.6s.mnn[logupx.6s.mnn < 0.05] <- 0
object.6s.mnn@logupx.data <- as(logupx.6s.mnn, "dgCMatrix")
#re-calculate PCA
object.6s.mnn <- calcPCA(object.6s.mnn)
#re-calculate tSNE
set.seed(18)
object.6s.mnn<- calcTsne(object.6s.mnn, dim.use= "pca", perplexity=30, theta=0.5)
#look at batch information
plotDim(object.6s.mnn, "BATCH", plot.title = "BATCH (MNN Corrected)")
#do graph clustering with Louvain-Jaccard
object.6s.mnn <- graphClustering(object.6s.mnn, num.nn = c(5,8,10,15), method = "Louvain", do.jaccard = T)
#do graph clustering with Infomap-Jaccard
object.6s.mnn <- graphClustering(object.6s.mnn, num.nn = c(10,15,20,30,40), method = "Infomap", do.jaccard = T)
clusterings <- c(paste0("Louvain-", c(5,8,10,15)), paste0("Infomap-", c(10,15,20,30,40)))
for (c in clusterings){plot(plotDim(object.6s.mnn, c, legend=F))}

clusters <- unique(object.6s.mnn@group.ids$`Infomap-30`)
#precision-recall markers to find the best 'markers' of each cluster
pr.markers <- lapply(clusters,function(c) markersAUCPR(object.6s.mnn, clust.1=c, clustering= "Infomap-30", genes.use= object.6s.mnn@var.genes))
names(pr.markers) <- clusters

#make a set of data.frames to keep track during cluster assignment
I30.n <- length(unique(object.6s.mnn@group.ids$`Infomap-30`))
I30.cluster.assignments <- data.frame(cluster=1:I30.n, name =rep(NA, I30.n), tip=rep(NA, I30.n), row.names= 1:I30.n)
I20.n <- length(unique(object.6s.mnn@group.ids$`Infomap-20`))
I20.cluster.assignments <- data.frame(cluster=1:I20.n, name =rep(NA, I20.n), tip=rep(NA, I20.n), row.names= 1:I20.n)

#--------------------------------------------------
#endoderm
plotDot(object.6s.mnn, genes =c("ENSDARG00000036292", "ENSDARG00000002601", "ENSDARG00000021232", "ENSDARG00000003411", "ENSDARG00000102138", "ENSDARG00000055064"), clustering = "Infomap-30")
I30.cluster.assignments["21","name"] <- "Endoderm Pharyngeal"
I30.cluster.assignments["40","name"] <- "Endoderm Pancreatic/Intestinal"

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "40", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000036292"), plotDim(object.6s.mnn, "ENSDARG00000055064"), plotDim(object.6s.mnn, "ENSDARG00000003411", legend=F)))

#axial mesoderm
plotDot(object.6s.mnn, genes =c("ENSDARG00000040944", "ENSDARG00000021201", "ENSDARG00000101576", "ENSDARG00000077403", "ENSDARG00000069093", "ENSDARG00000071082", "ENSDARG00000019122", "ENSDARG00000023656", "ENSDARG00000039173"), clustering = "Infomap-30")
I30.cluster.assignments["34","name"] <- "Prechordal Plate"
I30.cluster.assignments["43","name"] <- "Notochord Anterior" 

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "43", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000040944"), 
                        plotDim(object.6s.mnn, "ENSDARG00000101576"), 
                        plotDim(object.6s.mnn, "ENSDARG00000071082", legend=F)))


#intermediate/lateral mesoderm
plotDot(object.6s.mnn, genes =c("ENSDARG00000008305", "ENSDARG00000101919", "ENSDARG00000028148", "ENSDARG00000000767", "ENSDARG00000043271", "ENSDARG00000013477", "ENSDARG00000095019", "ENSDARG00000019930", "ENSDARG00000017195", "ENSDARG00000052846"), clustering = "Infomap-30")

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-20", "52", legend=F), plotDim(object.6s.mnn, "ENSDARG00000013477"), plotDim(object.6s.mnn, "ENSDARG00000043271"), plotDimHighlight(object.6s.mnn, "Infomap-10", "85", legend=F)))

plotDot(object.6s.mnn, genes =c("ENSDARG00000008305", "ENSDARG00000101919", "ENSDARG00000028148", 
                                "ENSDARG00000000767", "ENSDARG00000043271", "ENSDARG00000013477", 
                                "ENSDARG00000095019", "ENSDARG00000019930", "ENSDARG00000017195", 
                                "ENSDARG00000052846"), clustering = "Infomap-20")

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "7", legend=F),
plotDimHighlight(object.6s.mnn, "Infomap-20", "52", legend=F),
plotDimHighlight(object.6s.mnn, "Infomap-20", "53", legend=F),
plotDimHighlight(object.6s.mnn, "Infomap-20", "54", legend=F)
))
grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "7", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000008305"), 
                        plotDim(object.6s.mnn, "ENSDARG00000028148"), 
                        plotDim(object.6s.mnn, "ENSDARG00000000767", legend=F)))


I30.cluster.assignments["18","name"] <- "Cephalic Mesoderm"
I30.cluster.assignments["49","name"] <- "Hematopoietic (RBI)"
I20.cluster.assignments["52","name"] <- "Hematopoietic (ICM)"
I30.cluster.assignments["39","name"] <- "Pronephros"
I30.cluster.assignments["7","name"] <- "Heart Primordium"

#Paraxial mesoderm, no wnt8a
plotDot(object.6s.mnn, genes =c("ENSDARG00000070913", "ENSDARG00000101576", "ENSDARG00000021201", 
                                "ENSDARG00000003399", "ENSDARG00000006939","ENSDARG00000011785", 
                                "ENSDARG00000030347", "ENSDARG00000070535","ENSDARG00000053493",
                                "ENSDARG00000054103","ENSDARG00000007891","ENSDARG00000007277",
                                "ENSDARG00000009438","ENSDARG00000030110","ENSDARG00000062592"), clustering = "Infomap-30")

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "16", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000007277"), 
                        plotDim(object.6s.mnn, "ENSDARG00000007891"), 
                        plotDim(object.6s.mnn, "ENSDARG00000054103", legend=F)))


I30.cluster.assignments["48","name"] <- "Adaxial Cells"
#markers.11v33 <- markersAUCPR(object.6s.mnn, clust.1= "11", clust.2= "33", clustering = "Infomap-30")
I30.cluster.assignments["16","name"] <-"Somites Formed"
I30.cluster.assignments["5","name"] <- "Tailbud"
I30.cluster.assignments["37","name"] <- "Tailbud"

#Neural
plotDot(object.6s.mnn, genes =c("ENSDARG00000003411", "ENSDARG00000038867", "ENSDARG00000068567", 
                                "ENSDARG00000101919", "ENSDARG00000021032","ENSDARG00000043923", 
                                "ENSDARG00000077467"), clustering = "Infomap-30")
#no ENSDARG00000021032 Foxd3
I30.cluster.assignments["12","name"] <- "Neural crest"
I30.cluster.assignments["44","name"] <- "Floor Plate"

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "44", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000003411"), 
                        plotDim(object.6s.mnn, "ENSDARG00000038867"), 
                        plotDim(object.6s.mnn, "ENSDARG00000101919", legend=F)))

#spinal cord
plotDot(object.6s.mnn, genes =c("ENSDARG00000006350", "ENSDARG00000025017", "ENSDARG00000052610", 
                                "ENSDARG00000010791", "ENSDARG00000056130", "ENSDARG00000003469", 
                                "ENSDARG00000019566","ENSDARG00000053499","ENSDARG00000014420"), clustering = "Infomap-30")

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "28", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000010791"), 
                        plotDim(object.6s.mnn, "ENSDARG00000056130"), 
                        plotDim(object.6s.mnn, "ENSDARG00000003469", legend=F)))



I30.cluster.assignments["52","name"] <- "Spnial Cord Differentiated"
I30.cluster.assignments["28","name"] <- "Spinal Cord"
#I30.cluster.assignments["5","name"] <- "Spinal Cord Progenitors"

#Fore/mid-brain
plotDot(object.6s.mnn, genes =c("ENSDARG00000045936", "ENSDARG00000040321", "ENSDARG00000052893", 
                                "ENSDARG00000074253", "ENSDARG00000029179", "ENSDARG00000103379", 
                                "ENSDARG00000057936", "ENSDARG00000101549","ENSDARG00000001859",
                                "ENSDARG00000086393","ENSDARG00000020417","ENSDARG00000070769",
                                "ENSDARG00000028148","ENSDARG00000002707","ENSDARG00000008796",
                                "ENSDARG00000026599","ENSDARG00000038868"), clustering = "Infomap-30")
grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "10", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000045936"), 
                        plotDim(object.6s.mnn, "ENSDARG00000040321"), 
                        plotDim(object.6s.mnn, "ENSDARG00000052893", legend=F)))


I30.cluster.assignments["8","name"] <- "Midbrain"
I30.cluster.assignments["19","name"] <- "Telencephalon"
I30.cluster.assignments["35","name"] <- "Dorsal Diencephalon"
#I30.cluster.assignments["7","name"] <- "Dorsal Diencephalon"
I30.cluster.assignments["10","name"] <- "Optic cup"
I30.cluster.assignments["17","name"] <- "Optic cup"

#Hindbrain
plotDot(object.6s.mnn, genes =c("ENSDARG00000020164", "ENSDARG00000004782", "ENSDARG00000003399", 
                                "ENSDARG00000017121", "ENSDARG00000042826", "ENSDARG00000059276", 
                                "ENSDARG00000029263","ENSDARG00000103862","ENSDARG00000000175",
                                "ENSDARG00000023031","ENSDARG00000008174"), clustering = "Infomap-20")

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-20", "27", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000059276"), 
                        plotDim(object.6s.mnn, "ENSDARG00000029263"), 
                        plotDim(object.6s.mnn, "ENSDARG00000103862", legend=F)))


I20.cluster.assignments["28","name"] <- "Hindbrain R5+6"
I20.cluster.assignments["4","name"] <- "Hindbrain R3"
I20.cluster.assignments["35","name"] <- "Hindbrain R4"
I20.cluster.assignments["22","name"] <- "Hindbrain R4"
I20.cluster.assignments["27","name"] <- "Hindbrain R7"
I20.cluster.assignments["11","name"] <- "Hindbrain R7"

#Non-neural ectoderm
plotDot(object.6s.mnn, genes =c("ENSDARG00000102981", "ENSDARG00000006120", "ENSDARG00000059327", 
                                "ENSDARG00000073978", "ENSDARG00000014626", "ENSDARG00000045413", 
                                "ENSDARG00000053666","ENSDARG00000009550","ENSDARG00000055926"), clustering = "Infomap-30")

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "1", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000102981"), 
                        plotDim(object.6s.mnn, "ENSDARG00000006120"), 
                        plotDim(object.6s.mnn, "ENSDARG00000014626", legend=F)))


I30.cluster.assignments["50","name"] <- "Integument"
I30.cluster.assignments["15","name"] <- "Neural Plate Border"
I30.cluster.assignments["1","name"] <- "Epidermis"

#Pre-placodal ectoderm
plotDot(object.6s.mnn, genes =c("ENSDARG00000038792", "ENSDARG00000101831", "ENSDARG00000010477", 
                                "ENSDARG00000028148", "ENSDARG00000076480", "ENSDARG00000006120",
                                "ENSDARG00000104566","ENSDARG00000015879","ENSDARG00000060734",
                                "ENSDARG00000071560","ENSDARG00000054879","ENSDARG00000042525",
                                "ENSDARG00000010591","ENSDARG00000062892","ENSDARG00000116608",
                                "ENSDARG00000070069"), clustering = "Infomap-20")
#ENSDARG00000099564, ATOH1B

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-20", "50", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000038792"), 
                        plotDim(object.6s.mnn, "ENSDARG00000101831"), 
                        plotDim(object.6s.mnn, "ENSDARG00000010477", legend=F)))





I20.cluster.assignments["38","name"] <- "Placode Lens"
I20.cluster.assignments["79","name"] <- "Placode Olfactory"
I20.cluster.assignments["47","name"] <- "Placode Adenohypophyseal" 
I20.cluster.assignments["60","name"] <- "placode Epibranchial"
I20.cluster.assignments["51","name"] <- "Placode Otic"
I20.cluster.assignments["50","name"] <- "Placode Trigminal"

#NON-BLASTODERM
plotDot(object.6s.mnn, genes =c("ENSDARG00000068255",  "ENSDARG00000022813", "ENSDARG00000014373", 
                                "ENSDARG00000093864", "ENSDARG00000018404", "ENSDARG00000036834",
                                "ENSDARG00000017624","ENSDARG00000058371","ENSDARG00000094041"), clustering = "Infomap-30")
#ENSDARG00000102892, H1M IS NOT ANNOTATED

grid.arrange(grobs=list(plotDimHighlight(object.6s.mnn, "Infomap-30", "55", legend=F), 
                        plotDim(object.6s.mnn, "ENSDARG00000068255"), 
                        plotDim(object.6s.mnn, "ENSDARG00000022813"), 
                        plotDim(object.6s.mnn, "ENSDARG00000014373", legend=F)))


evl.score<-apply(object@logupx.data[c("ENSDARG00000093864", "ENSDARG00000018404", "ENSDARG00000036834",
                                      "ENSDARG00000017624","ENSDARG00000058371","ENSDARG00000094041"), 
                                    cellsInCluster(object.6s.mnn, "Infomap-30", "38")],2, sum.of.logs)

new.evl <-names(which(evl.score >9))
remove.evl <-names(which(evl.score <=9))
I30.cluster.assignments["38","name"] <- "EVL/Periderm"
I30.cluster.assignments["55","name"] <- "Primordial Germ Cells"
#I30.cluster.assignments["64","name"] <- "YSL"

#general final clustering
I30.cluster.assignments$clustering <-"Infomap-30"
I20.cluster.assignments$clustering <-"Infomap-20"
cluster.assignments <- rbind(I30.cluster.assignments, I20.cluster.assignments)
cluster.assignments <- cluster.assignments[!is.na(cluster.assignments$name), ]
cluster.assignments$cluster.new <- 1:nrow(cluster.assignments)
object.6s.mnn@group.ids$clusters.6s.name <-NA
object.6s.mnn@group.ids$clusters.6s.num <-NA

#copy cell identities over for each cluster  (not sure about the for loop)
for (i in 1:nrow(cluster.assignments)){
  cells<- cellsInCluster(object.6s.mnn, clustering = cluster.assignments[i, "clustering"], cluster = cluster.assignments[i, "cluster"]) 
  object.6s.mnn@group.ids[cells, "clusters.6s.name"] <- cluster.assignments[i, "name"] 
  object.6s.mnn@group.ids[cells, "clusters.6s.num"] <-as.character(cluster.assignments[i, "cluster.new"])
}
clusters

object.6s.mnn@group.ids[remove.evl, "cluster.6s.name"] <- NA
object.6s.mnn@group.ids[remove.evl, "cluster.6s.num"] <- NA

object@group.ids$`ZF6S-Infomap-30` <- NA
object@group.ids[rownames(object.6s.mnn@group.ids), "ZF6S-Infomap-30"] <- object.6s.mnn@group.ids$`Infomap-30`

object@group.ids$`ZF6S-Infomap-20` <- NA
object@group.ids[rownames(object.6s.mnn@group.ids), "ZF6S-Infomap-20"] <-object.6s.mnn@group.ids$`Infomap-20`

object@group.ids$`ZF6S-Cluster` <- NA
object@group.ids[rownames(object.6s.mnn@group.ids), "ZF6S-Cluster"] <-object.6s.mnn@group.ids$clusters.6s.name

 FALSE)
  }
}
############################################################
object<-readRDS("object_4_withTips_dm11_0910.rds")
diffusion.logistic <- pseudotimeDetermineLogistic(object,
                                                  "pseudotime",
                                                  optimal.cells.forward = 40,
                                                  max.cells.back = 100,
                                                  pseudotime.direction = "<",
                                                  do.plot = TRUE,
                                                  print.values = TRUE)

# Bunch of stuff gets done on the cluster here.

tip.walk.files <- list.files("./walks_40_100/dm-11-tm-40-100/ZF6S-Cluster-Num/", pattern= ".rds",
                             full.names = TRUE)
tips.walked <- setdiff(unique(object@group.ids$`ZF6S-Cluster-Num`), NA)
for (tip in tips.walked) {
  tip.files <- grep(paste0("walks-", tip, "-"), tip.walk.files, value = TRUE)
  if (length(tip.files) > 0) {
    these.walks <- unlist(lapply(tip.files, readRDS), recursive = FALSE)
    object <<- processRandomWalks(object,
                                  walks = these.walks,
                                  walks.name = tip,
                                  verbose = FALSE)
  }
}

saveRDS(object, file= "object_5_with_walks_40_100.rds")

#############################################################
object <- readRDS("object_5_with_walks_40_100.rds")
object <- loadTipCells(object, tips = "ZF6S-Cluster-Num")
#tip.clusters <- read.csv("tips-use_6s_dm11_0910.csv")

#hindbrain R4
object<- combineTipVisitation(object, "29", "32", "29")
object<- combineTipVisitation(object, "28", "30","28")
object<- combineTipVisitation(object, "5", "9","5")
object<- combineTipVisitation(object, "2", "16","16")
#endo
object<- combineTipVisitation(object, "12", "19","12")

#hermatopoietic
object<- combineTipVisitation(object, "23", "37","23")

tips.to.exclude<-c("32","30","9","2","19","37")
#didnt exclude tips that didnt use
tips.to.use <-setdiff(as.character(1:39), tips.to.exclude)

tips.to.use <-as.character(object@tree$tips)
tips.to.use


object.built <- buildTree(object = object,
                          pseudotime = "pseudotime",
                          divergence.method = "ks",
                          tips.use = tips.to.use,
                          weighted.fusion = TRUE,
                          use.only.original.tips = TRUE,
                          cells.per.pseudotime.bin = 80,
                          bins.per.pseudotime.window = 6,
                          minimum.visits = 1,
                          visit.threshold = 0.7,
                          p.thresh = 0.025,
                          save.breakpoint.plots = NULL,
                          dendro.node.size = 100,
                          min.cells.per.segment = 10,
                          min.pseudotime.per.segment = 0.01,
                          verbose = TRUE)

tip.names <- unique(object@group.ids[, c("ZF6S-Cluster", "ZF6S-Cluster-Num")])
tip.names <- tip.names[complete.cases(tip.names), ]
object.built <- nameSegments(object.built, 
                             segments = tip.names$`ZF6S-Cluster-Num`,
                             segment.names = tip.names$`ZF6S-Cluster`)
plotTree(object.built, label.segments=T)
#test for gene/tes
plotTree(object.built, "ERV2-DR-LTR")
plotTree(object.built, "KibiDr1")
plotTree(object.built, "BHIKHARI-LTR")
plotTree(object.built, "ENSDARG00000020417")

#change names
#new.seg.names <- c("Optic cup", "Heart Primordium", "Midbrain", "Neural crest", "Somites Formed", "Dorsal Diencephalon", "Epidermis", "Neural Plate Border", "Telencephalon", "Cephalic Mesoderm", "Endoderm", "Hematopoietic", "Prechordal Plate", "Pronephros", "EVL/Periderm", "Notochord Anterior", "Adaxial Cells", "Integument", "Floor Plate", "Spnial Cord Differentiated", "Primordial Germ Cells", "Placode Trigminal", "Placode Lens", "Placode Otic", "Placode Adenohypophyseal", "Placode Olfactory", "Hindbrain")
#new.short.names <- c("Optic","Heart","Mid", "NC","Somite","Dien","Epi","NPB","Tel","CM","Endo","Hem","PP","Pro","EVL","Noto","Adax","Int","FP","SC","PGC","PT","PL","POtic","PA","POlf","Hind")
#segs.to.name <- c("1","2","3","4","5","6","7","8","9","10","12","13","15","17","18","19","21","22","23","24","25","27","30","31","32","34","35")
new.seg.names <- c("Optic cup","Heart Primordium","Midbrain","Neural crest","Somites Formed","Dorsal Diencephalon","Epidermis","Neural Plate Border","Cephalic Mesoderm","Endoderm","Hematopoietic","Prechordal Plate","Pronephros","EVL/Periderm","Notochord Anterior","Adaxial Cells","Integument+Placode","Floor Plate","Spnial Cord Differentiated","Primordial Germ Cells","Hindbrain R7","Hindbrain R5+6","Hindbrain R4","Tailbud")
new.short.names <- c("Optic","Heart","Mid","NC","Somite","Dien","Epi","NPB","CM","Endo","Hem","PP","Pro","EVL","Noto","Adax","Int+PC","FP","SC","PGC","Hind7","Hind5+6","Hind4","TB")
segs.to.name <- c("5","3","47","42","8","48","1","7","10","12","23","14","18","17","20","22","44","45","41","26","28","46","29","16")

object.built<-nameSegments(object.built, segments=segs.to.name, segment.names=new.seg.names, short.names=new.short.names)
plotTree(object.built, label.segments=T)



#create force layout
#Data frame to measure cell visitation
visitation <-data.frame(cell=rownames(object.built@diff.data), 
                        seg=object.built@diff.data$segment,
                        stringsAsFactors=F,
                        row.names=rownames(object.built@diff.data))
visitation$visit<-log10(apply(visitation, 1, function(cr) object.built@diff.data[as.character(cr["cell"]), paste0("visitfreq.raw.", as.character(cr["seg"]))]) +1)

#Choose those cells that were well visited
robustly.visited.cells <- visitation[visitation$visit >=0.5, "cell"]

final.tips <- segTerminal(object.built)
#generate the force-directed layout
object.built <- treeForceDirectedLayout(object.built, num.nn=120, 
                                          method= "fr", 
                                          cells.to.do =robustly.visited.cells, 
                                          tips =final.tips, 
                                          cut.unconnected.segments=2, 
                                          min.final.neighbors =4, verbose =T)

#no product on the server, save it anyway
saveRDS(object, file= "object_5.5_half_tree_40_100.rds")
saveRDS(object.built, file= "object_5.5_object.built.rds")



#load view
#object.built@tree$force.view.list <- readRDS("force.view.list.rds")
#object.built@tree$force.view.default <- "figure1"
#load layout
#precalc.fdl <-readRDS("force.view.list.rds")
#object.built@tree$walks.force.layout <- precalc.fdl$walks.force.layout
plotTreeForce(object.built.2, "BHIKHARI-LTR",alpha =0.2, view="figure1") #view= "figure1"

#hand-tune the tree
object.built.3 <- treeForceRotateCoords(object.built.2, seg ="52", angle= -96, axis = "z", around.cell =10, throw.out.cells =1000, pseudotime= "pseudotime")

#curve
for (throw.out in c(0, 100, 250, 500, 1000)) {
  object.built.3 <- treeForceRotateCoords(object.built.3, seg= "52", angle = pi/5, axis="y", around.cell =10, throw.out.cells = throw.out, pseudotime= "pseudotime")
}

for (throw.out in c(0, 100, 250, 500, 1000)) {
  object.built.3 <- treeForceRotateCoords(object.built.3, seg= "52", angle = pi/10, axis="z", around.cell =10, throw.out.cells = throw.out, pseudotime= "pseudotime")
}

plotTreeForce(object.built.3, "BHIKHARI-LTR",alpha =0.2, view="figure1")
#axial mesoderm
object.built.4 <- treeForceRotateCoords(object.built.3, seg ="56", angle= -36, axis = "y", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.4 <- treeForceRotateCoords(object.built.4, seg ="56", angle= 36, axis = "x", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.4 <- treeForceRotateCoords(object.built.4, seg ="56", angle= 20, axis = "z", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")

plotTreeForce(object.built.4, "BHIKHARI-LTR",alpha =0.2, view="figure1")

object.built.5 <- treeForceRotateCoords(object.built.5, seg ="55", angle= 3, axis = "y", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.5 <- treeForceRotateCoords(object.built.4, seg ="55", angle= 5, axis = "x", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.5 <- treeForceRotateCoords(object.built.5, seg ="55", angle= 3, axis = "z", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")

plotTreeForce(object.built.5, "BHIKHARI-LTR",alpha =0.2, view="figure1")


#remainder of the mesendoderm
object.built.6 <- treeForceRotateCoords(object.built.6, seg ="54", angle= -12, axis = "z", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.6 <- treeForceRotateCoords(object.built.6, seg ="54", angle= -5, axis = "y", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.6 <- treeForceRotateCoords(object.built.6, seg ="54", angle= -5, axis = "x", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")

plotTreeForce(object.built.6, "BHIKHARI-LTR",alpha =0.2, view="figure1")

#try to bend just one branch
#heart
object.built.7 <- treeForceRotateCoords(object.built.7, seg ="2", angle= 13, axis = "z", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.7 <- treeForceRotateCoords(object.built.7, seg ="2", angle= -13, axis = "x", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.7 <- treeForceRotateCoords(object.built.7, seg ="2", angle= -6, axis = "y", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")

plotTreeForce(object.built.7, "BHIKHARI-LTR",alpha =0.2, view="figure1")

#CM 10
object.built.8 <- treeForceRotateCoords(object.built.8, seg ="10", angle= -28, axis = "z", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.8 <- treeForceRotateCoords(object.built.8, seg ="10", angle= -10, axis = "x", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.8 <- treeForceRotateCoords(object.built.8, seg ="10", angle= -3, axis = "y", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")

plotTreeForce(object.built.8, "BHIKHARI-LTR",alpha =0.2, view="figure1")

#hind
object.built.9 <- treeForceRotateCoords(object.built.8, seg ="35", angle= -28, axis = "z", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.9 <- treeForceRotateCoords(object.built.8, seg ="35", angle= 0, axis = "x", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.9 <- treeForceRotateCoords(object.built.9, seg ="35", angle= 3, axis = "y", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")

plotTreeForce(object.built.9, "BHIKHARI-LTR",alpha =0.2, view="figure1")



#germ
object.built.10 <- treeForceRotateCoords(object.built.9, seg ="25", angle= -pi/2, axis = "y", around.cell =1, throw.out.cells =0, pseudotime= "pseudotime")
object.built.10 <- treeForceRotateCoords(object.built.10, seg ="25", angle= 183, axis = "z", around.cell =1, throw.out.cells =0, pseudotime= "pseudotime")
object.built.10 <- treeForceRotateCoords(object.built.10, seg ="25", angle= 5, axis = "x", around.cell =1, throw.out.cells =0, pseudotime= "pseudotime")


object.built.12 <- treeForceTranslateCoords(object.built.11, seg= "25", x = 0, y = 0, z = -3)

plotTreeForce(object.built.12, "BHIKHARI-LTR",alpha =0.2, view="figure1")

#re-tweet hind, otherwise is fine
object.built.11 <- treeForceRotateCoords(object.built.10, seg ="35", angle= -3, axis = "z", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.11 <- treeForceRotateCoords(object.built.11, seg ="35", angle= 0, axis = "x", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")
object.built.11 <- treeForceRotateCoords(object.built.11, seg ="35", angle= 3, axis = "y", around.cell =10, throw.out.cells =0, pseudotime= "pseudotime")


library(rgl)
library(gridExtra)
library(RColorBrewer)
rgl::setupKnitr()
object12 <-readRDS("object.built.12.rds")
branch.colors <- c("#CECECE","#E6298B")
fire.with.grey <- c("#CECECE","#DDC998",RColorBrewer::brewer.pal(9,"YlOrRd")[3:9])
plotTreeForce(object12, "ENSDARG00000052846",alpha =0.5,
              colors = fire.with.grey, label.tips=FALSE)
ggsave(plot, filename="test.pdf")
EnSpmN1-DR

#Determine EVL cells
#evl.cells <-intersect(cellsInCluster(object.built, "segment","38"), rownames(object.built@tree$walks.force.layout))
#evlcells.move.1
#evl.cells.move.2
#object.built <- treeForceTranslateCoords(object.built, cells=evl.cells.move.1, x=0,y=0,z=-10)
#object.built <- treeForceTranslateCoords(object.built, cells=evl.cells.move.2, x=0,y=0,z=-15)

object.built.3 <- treeForceRotateCoords(object.built.3, seg ="18", angle= 1.35, axis = "y", around.cell =1, throw.out.cells =0, pseudotime= "pseudotime")

plotTreeForce(object.built.3, "segment", alpha= 0.2, view ="figure1")
```
