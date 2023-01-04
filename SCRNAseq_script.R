library(ggplot2)
#BiocManager::install('Seurat')
library(Seurat)
library(dplyr)
#BiocManager::install("celda") 
library(celda)
#library(celda) the updataed celda has all new functions 
#install.packages("usethis")
#library(usethis)
#install.packages("devtools")
#library(devtools)
#install_github("campbio/celda@20190409_master")
data = read.csv("/home/genomics/single_cell_data/RNA_umi.csv", header = TRUE,  row.names = 1)

data_protein = read.csv("/home/genomics/single_cell_data/Protein_umi.csv", header = TRUE,  row.names = 1)
#head(data, n=10) 
#create an assay object for multimodal analysis
alldata <- CreateSeuratObject(counts = data)                
alldata[["protein"]] <-CreateAssayObject(counts = data_protein)
Assays(alldata)
alldata <- NormalizeData(alldata, assay = "RNA", normalization.method = "CLR")
#this counts the number of cells that are not empty in each column (cell)
#numberofexpressed = as.data.frame(colSums(data!=0))
#we want to show he number of genes expressed in each cell, this will help us understand the quality
#ggplot(numberofexpressed) +
#  geom_histogram(
#    breaks = seq(1000,10000, by = 100),
#    aes(x=(colSums(data!=0))))
#scale_x_continuous(breaks = seq(1000,10000, by = 100))
data_test <- data[1:36280, 1:10]
#100 cells takes 2mins, 8000 cells will take 5hrs! 2seconds per cell
library(Matrix)
mat2 <- Matrix(mat)
mat <- as.matrix(data_test)
mat_data <- as.matrix(human_clusters@"assay"[["RNA"]]@"data")

#  $human_clusters@"assays"[["RNA"]]@"data")
storage.mode(mat) <- "logical"
storage.mode(mat_data) <-"logical"
#BiocManager::install("SingleCellExperiment")
#library(SingleCellExperiment)

#data_sce <- SingleCellExperiment(mat)

data_decont <-decontX(mat)

#can export the top30 genes in each cell type and run pathway enrichment analysis in Gene ontology
#DEseq can help to show differential gene expression
BiocManager::install("DESeq2") 
#we will need to get a coldata.csv from information on cells in the count matrix
#rows in coldatwe should correspind to cols in the count matrix
#maybe we can get from the cellTypeMappings
colData <- colData(data_decontx)
dds <- DESeqDataSetFromMatrix(countData = data_decontx,
                              colData = coldata(data_decontx),
                              design = ~ cell + type)

nrows(dds)

keep<-rowSums(counts(dds)) >20
dds<- dds[keep,]
nrows
vsd <- vst(dds, blind = FALSE)
#runs the deseq and plots results
sampleDist <- dist(t(assay(vsd)))
dds <-DESeq(dds)
res <- results(dds)

sum(res$padj <0.1, na.rm = TRUE)
resSig <-subset(res,padj <0.1)

head(resSig[order(resSig$log2FoldChange, decreasing = TRUE),])

library("apeglm")
resutsNames(dds)

res <- lfcShrink(dds, coef "Monocytes")
#cellRouters can help to show cell fata trajectories

#save.image(file = "scRNAseq_PBMC")
#creates an object of Seurat class from count matrix from library Seurat

alldata, assay = 'RNA", <- CreateSeuratObject(counts = data, project = "PBMC", min.cells = 3, min.features= 2)

# a lot of documentation suggests to put "^MT-" it didn't work here, "MT-" works.
#grep("MT-", rownames(data_Seurat@assays$RNA@counts), value = TRUE)
#adds the percent of mitochondiral genes in each cell, a live sample will have less than 15%          
data_Seurat[["percent.mt"]] <- PercentageFeatureSet(object = data_Seurat, pattern = "MT-")
data_Seurat[["percent.mouse"]] <- PercentageFeatureSet(object = data_Seurat, pattern = "MOUSE")

#visualise QC metrics following checks makes sure the data is in the suerat file 
#and that the datis in the right formats 

#grep("^MT", rownames(data_Seurat[["RNA"]]), value = T)

#data_Seurat@assays$RNA@counts[1:10,1:10]

VlnPlot(data_Seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.2)
VlnPlot(data_Seurat, features=c("percent.mt","percent.mouse"), ncol = 2, pt.size = 0.2)


FeatureScatter(data_Seurat, feature1= "nCount_RNA", feature2="percent.mt")
FeatureScatter(data_Seurat, feature1= "nCount_RNA", feature2="nFeature_RNA")


#We need to normalise to make it easier to visualise the expression differences
data_Seurat <-NormalizeData(data_Seurat,normalization.method = "LogNormalize")

#We want to find the most variable genes 
data_Seurat <-FindVariableFeatures(data_Seurat, selection.method = "vst", nfeatures = 2000)
#topten <-head(VariableFeatures(data_Seurat),10)
#topten
#plot1 <-VariableFeaturePlot(data_Seurat)
#plot2 <- LabelPoints(plot= plot1, points = topten, repel = TRUE)
#plot2

all.genes <-rownames(data_Seurat)
data_Seurat <- ScaleData(data_Seurat, features = all.genes)
#We want to run a PCA on the data and show a few plots 

data_Seurat<-RunPCA(object = data_Seurat, features = VariableFeatures(object = data_Seurat))                         
#print(data_Seurat[["pca"]], dims = 1:2, nfeatures = 5)      

#VizDimLoadings(data_Seurat, dims = 1:2, nfeatures = 15, reduction = "pca")

DimPlot(data_Seurat, reduction ="pca")
#DimHeatmap(data_Seurat, dims = 1, cells =500, balanced = TRUE)

#we want to determine how many clusters we have in the data
ElbowPlot(data_Seurat)

data_Seurat<-FindNeighbors(data_Seurat, dims = 1:25)
data_Seurat<-FindClusters(data_Seurat, resolution = 1.04)

head(Idents(data_Seurat), 25)

#UMPA is better for visualising a lot of dat
data_Seurat <- RunUMAP(data_Seurat, dims = 1:25)
#when you run the below, it shouls have a graph with cells coloured by cluster
#if it doesn't  forgot to run find clusters
DimPlot(data_Seurat, reduction = "umap")

#Now we wnat to find what the clusters are, and which ones represent the mouse cells, later we can remove these clusters

#the following clusters are mouse: 0 ,1, 5, 8 need to remove and rescale
#can subset 
#VlnPlot(data_Seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.2)
#VlnPlot(data_Seurat, features= c("percent.mt", "percent.mouse"), ncol = 2, pt.size = 0.2)
# plots show that cluster 5 and 8 and really moue clusters, 1 and 7 have a few but we can clean using thresholds
# clean up the cells that have more than 1% mouse and <10% mt, and gene exppression of less than 50 genes
#we are going to cluster the gene expression that have top genes as mouse so we can filture those clustsers
data_Seurat_markers <- FindAllMarkers(data_Seurat, only.pos = TRUE, pin.pct = 0.25, logfc.threshold = 0.25)
#it takes aroun 10mins

top10 <- data_Seurat_markers %>%
  group_by(cluster)%>%
  top_n(n = 10)
top10
#this plots a heatmap with the top ten genes per cluster, we can make sure out clusters are good
p2 <- DoHeatmap(data_Seurat, features = top10$gene, group.bar.height = 0.01, size = 3, combine = FALSE)
p2 <-lapply(X=p2, FUN = function(x)x+
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5))+
              theme(axis.title.x = element_text(size = 5))+
              theme(axis.text.y = element_text(size = 3))+
              theme(axis.text.x = element_text(size = 3))+
              theme(legend.position = "none"))
CombinePlots(plots = p2)

# Note, for simplicity we are merging two CD14+ Monocyte clusters (that differ in expression of
# HLA-DR genes) and NK clusters (that differ in cell cycle stage)
new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B", 
                     "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", 
                     "Mouse", "DC", "pDCs")
names(new.cluster.ids) <- levels(data_Seurat)
data_Seurat <- RenameIdents(data_Seurat, new.cluster.ids)
DimPlot(data_Seurat, label = TRUE) + NoLegend()


data_Seurat[["protein"]] <-CreateAssayObject(counts = data_protein)
data_Seurat <- NormalizeData(data_Seurat, assay = "protein", normalization.method = "CLR")
#data_Seurat <- FindVariableFeatures(alldata, assay = "protein")
data_Seurat <- ScaleData(data_Seurat, assay = "protein")


p <- FeaturePlot(data_Seurat, features = c(
p <-lapply(X=p, FUN = function(x)x+
             theme(plot.title = element_text(size = 8)) +
             theme(axis.title.y = element_text(size = 5))+
             theme(axis.title.x = element_text(size = 5))+
             theme(axis.text.y = element_text(size = 5))+
             theme(axis.text.x = element_text(size = 5))+
             theme(legend.position = "none"))
CombinePlots(plots = p)                                      

human_clusters <- subset(x=data_Seurat, idents = c(ident.1 = 5, ident.1 = 8, ident.1 = 10, indent.1 =11), invert = TRUE)
human_clusters <- subset(x=human_clusters, subset = percent.mouse < 5)
human_clusters <- subset(x= human_clusters, subset  = percent.mt <10)
human_clusters <- subset(x = human_clusters, subset = nFeature_RNA > 50)
human_clusters[["percent.mt_cleaned"]] <- PercentageFeatureSet(object = human_clusters, pattern = "MT-")
human_clusters[["percent.mouse_cleaned"]] <- PercentageFeatureSet(object = human_clusters, pattern = "MOUSE")

#plot the new graphs to check
VlnPlot(human_clusters, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.2)
VlnPlot(human_clusters, features = c("percent.mt","percent.mouse"), ncol = 2, pt.size = 0.2)

#renormalise the data 
human_clusters <-FindVariableFeatures(human_clusters, selection.method = "vst", nfeatures = 2000)
topten <-head(VariableFeatures(human_clusters),10)
#topten
plot1 <-VariableFeaturePlot(human_clusters)
plot2 <- LabelPoints(plot= plot1, points = topten, xnudge = 0, ynudge = 0, repel =TRUE)
plot

all.genes <-rownames(human_clusters)
human_clusters <- ScaleData(human_clusters, features = all.genes)
#rerun the PCA again to make sure the clusters are ok

human_clusters<-RunPCA(object = human_clusters, features = VariableFeatures(object = human_clusters))                         
#print(human_clusters[["pca"]], dims = 1:2, nfeatures = 5)     

VizDimLoadings(human_clusters, dims = 1:2, nfeatures = 15, reduction = "pca")

DimPlot(human_clusters, reduction ="pca")
#change the dims to 1,2,3 to get heatmaps for the respsective PCS
#DimHeatmap(human_clusters, dims = 3, cells =500, balanced = TRUE)

#we want to determine how many clusters we have in the data

ElbowPlot(human_clusters)

human_clusters<-FindNeighbors(human_clusters, dims = 1:10)
human_clusters<-FindClusters(human_clusters, resolution = 0.2)

head(Idents(human_clusters), 16)

#UMPA is better for visualising a lot of dat
human_clusters <- RunUMAP(human_clusters, dims = 1:10)



DimPlot(human_clusters, reduction = "umap")

#find markers again but just do find all this time 
human_clusters.markers <- FindAllMarkers(human_clusters, only.pos = TRUE, pin.pct = 0.25, logfc.threshold = 0.25)


x <-human_clusters.markers %>% group_by(cluster) %>% top_n(n = 1)
FeaturePlot(human_clusters, features = x$gene[1:4])
FeaturePlot(human_clusters, features = x$gene[5:9])
FeaturePlot(human_clusters, features = x$gene[10:11])

p <- FeaturePlot(human_clusters, features = c("HUMAN-CD3D",
                                              "HUMAN-CD3E",
                                              "HUMAN-CD4",
                                              "HUMAN-IL7R",
                                              "HUMAN-CD8A"),
                                              combine = FALSE)
p <- lapply(X=p, FUN = function(x)x+
             theme(plot.title = element_text(size = 8))+
             theme(axis.title.y = element_text(size = 5))+
             theme(axis.title.x = element_text(size = 5))+
             theme(axis.text.y = element_text(size = 5))+
             theme(axis.text.x = element_text(size = 5))+
             theme(legend.position = "none"))
CombinePlots(plots = p)
#clusters 9,0,5,8 show signs of tcells
                                              
p2 <- FeaturePlot(human_clusters, features = c("HUMAN-GNLY", #NK cells
                                             "HUMAN-NKG7", #NK cells
                                             "HUMAN-CD79A", #bcells
                                             "HUMAN-CD79B", #bcells
                                             "HUMAN-MS4A1"),
                                             combine = FALSE)
p2 <-lapply(X=p2, FUN = function(x)x+
             theme(plot.title = element_text(size = 8)) +
             theme(axis.title.y = element_text(size = 5))+
             theme(axis.title.x = element_text(size = 5))+
             theme(axis.text.y = element_text(size = 5))+
             theme(axis.text.x = element_text(size = 5))+
             theme(legend.position = "none"))
CombinePlots(plots = p2)
clustes 2,3,9
p3 <-FeaturePlot(human_clusters, features = c("HUMAN-S100A9",
                                             "HUMAN-FCGR3A", #fcGR3A+ MONOCYTES
                                             "HUMAN-CD14", #CD14 + MONOCYTES
                                             "HUMAN-LYZ",
                                             combine = FALSE))
#

- FeaturePlot(human_clusters, features = c("MS4A1","GNLY","CD3E","CD14", "FCER1A", "FCGR3A","LYZ","PPBP","CD8A"), combine = FALSE)
p <-lapply(X=p, FUN = function(x)x+
             theme(plot.title = element_text(size = 8)) +
             theme(axis.title.y = element_text(size = 5))+
             theme(axis.title.x = element_text(size = 5))+
             theme(axis.text.y = element_text(size = 5))+
             theme(axis.text.x = element_text(size = 5))+
             theme(legend.position = "none"))
CombinePlots(plots = p)
#tCELL(2), NKcells(1), bcells (3), monocyte(3)

top10 <- human_clusters.markers %>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_logFC)
top10
#this plots a heatmap with the top ten genes per cluster, we can make sure out clusters are good
p2 <- DoHeatmap(human_clusters, features = top10$gene, group.bar.height = 0.01, size = 3, combine = FALSE)
p2 <-lapply(X=p2, FUN = function(x)x+
             theme(plot.title = element_text(size = 6)) +
             theme(axis.title.y = element_text(size = 5))+
             theme(axis.title.x = element_text(size = 5))+
             theme(axis.text.y = element_text(size = 3))+
             theme(axis.text.x = element_text(size = 3))+
             theme(legend.position = "none"))
CombinePlots(plots = p2)

p5 <-FeaturePlot(data_Seurat, features = c("CD4", "CD45RA",
                                           "CD8", "CD3", "CCR7",
                                           "CD19", "CD16", "CD56","CD14"), combine =FALSE)
p5 <-lapply(X=p5, FUN = function(x)x+
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5))+
              theme(axis.title.x = element_text(size = 5))+
              theme(axis.text.y = element_text(size = 5))+
              theme(axis.text.x = element_text(size = 5)+
              theme(legend.position = "none"))
CombinePlots(plots = p5)


tcells <- subset(human_clusters, features = (ident.1 = 0),)
tcells <-FindVariableFeatures(tcells, selection.method = "vst", nfeatures = 2000)


nkcells <- subset(human_clusters, features = c(ident.1 = 2, ident.1 = 9))
nkcells <-FindVariableFeatures(nkcells, selection.method = "vst", nfeatures = 2000)


bcells <- subset(human_clusters, features = c(ident.1 = 4))
bcells <-FindVariableFeatures(tcells, selection.method = "vst", nfeatures = 2000)

monocytes <- subset(human_clusters, features = c(ident.1 = 5, ident.1 = 1))
moncytes <-FindVariableFeatures(monocytes, selection.method = "vst", nfeatures = 2000)



