library(ggplot2)
#BiocManager::install('Seurat')
  library(Seurat)
  library(dplyr)
#BiocManager::install("celda") 
library(celda)
BiocManager::install("SeuratWrappers")
#library(celda) the updataed celda has all anew functions 
#install.packages("usethis")
#library(usethis)
#install.packages("devtools")
#library(devtools)
#install_github("campbio/celda@20190409_master")
data = read.csv("/home/genomics/single_cell_data/RNA_umi.csv", header = TRUE,  row.names = 1)

data_protein = read.csv("/home/genomics/single_cell_data/Protein_umi.csv", header = TRUE,  row.names = 1)
#head(data, n=10) 
#create an assay object for multimodal analysis
data_Seurat<- CreateSeuratObject(counts = data)      
rm(data)

alldata[["protein"]] <-CreateAssayObject(counts = data_protein)

rm(data_protein)
Assays(data_Seurat)
data_Seurat <- NormalizeData(data_Seurat, assay = "RNA", normalization.method = "CLR")
data_Seurat <- FindVariableFeatures(data_Seurat, assay = "RNA")
data_Seurat <- ScaleData(data_Seurat, assay = "RNA")

data_Seurat <- NormalizeData(data_Seurat, assay = "protein", normalization.method = "CLR")
data_Seurat <- FindVariableFeatures(data_Seurat, assay = "protein")
data_Seurat<- ScaleData(data_Seurat, assay = "protein")
#this counts the number data_Seurat cells that are not empty in each column (cell)
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
colData <- colData(data_Seurat)
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

#data_Seurat<- CreateSeuratObject(counts = data, project = "PBMC", min.cells = 3, min.features= 2)

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
data_Seurat <-NormalizeData(data_Seurat, assay = "RNA", normalization.method = "LogNormalize")

#We want to find the most variable genes 
data_Seurat <-FindVariableFeatures(data_Seurat, assay = "RNA", selection.method = "vst", nfeatures = 2000)
#topten <-head(VariableFeatures(data_Seurat),10)
#topten
#plot1 <-VariableFeaturePlot(data_Seurat)
#plot2 <- LabelPoints(plot= plot1, points = topten, repel = TRUE)
#plot2

all.genes <-rownames(data_Seurat)
data_Seurat <- ScaleData(data_Seurat, assay = "RNA", features = all.genes)
#We want to run a PCA on the data and show a few plots 

data_Seurat<-RunPCA(object = data_Seurat, assay = "RNA", features = VariableFeatures(object = data_Seurat, assay = "RNA"))                         
#print(data_Seurat[["pca"]], dims = 1:2, nfeatures = 5)      

#VizDimLoadings(data_Seurat, dims = 1:2, nfeatures = 15, reduction = "pca")

DimPlot(data_Seurat, reduction ="pca")
#DimHeatmap(data_Seurat, )

#we want to determine how many clusters we have in the data
ElbowPlot(data_Seurat)

data_Seurat<-FindNeighbors(data_Seurat, assay = "RNA", dims = 1:25)
data_Seurat<-FindClusters(data_Seurat, assay = "RNA", resolution = 1.04)
#try to get 19 clusters
head(Idents(data_Seurat, assay = "RNA"), 25)

#UMPA is better for visualising a lot of dat
data_Seurat <- RunUMAP(data_Seurat, assay = "RNA", dims = 1:25)
#when you run the below, it shoFuls have a graph with cells coloured by cluster
#if it doesn't  forgot to run find clusters

DimPlot(data_Seurat)


#the following clusters are mouse: 0 ,1, 5, 8 need to remove and rescale
#can subset 
#VlnPlot(data_Seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.2)
#VlnPlot(data_Seurat, features= c("percent.mt", "percent.mouse"), ncol = 2, pt.size = 0.2)
# plots show that cluster 5 and 8 and really moue clusters, 1 and 7 have a few but we can clean using thresholds
# clean up the cells that have more than 1% mouse and <10% mt, and gene exppression of less than 50 genes
#we are going to cluster the gene expression that have top genes as mouse so we can filture those clustsers
data_Seurat_markers <- FindAllMarkers(data_Seurat, only.pos = TRUE, pin.pct = 0.25, logfc.threshold = 0.25)
#it takes aroun 10mins

# Note, for simplicity we are merging two CD14+ Monocyte clusters (that differ in expression of
# HLA-DR genes) and NK clusters (that differ in cell cycle state
#becausethe data set in UMAP orders the clusters by size, we can reuse someone elses idents for the cell types, 
#as long as we have the same number of clusurs. 

new.cluster.ids <- c("0", "1","2","3","4",
                    "5", "6", "7","8","9","10","11","12","13","14","15","16","17","18")
new.cluster.ids <- (c("Naive CD4+ T", "T/Mono doublets","CD14+ Mono", "Eryth","mouse","CD16+ Mono","Multiplets", "Memory CD4 T","NK", "CD8+ T", "pDCs","DCs", "Multiplets","Multiplets", "Multiplets","Multiplets"," B", "MK","Multiplets"))
names(new.cluster.ids) <- levels(data_Seurat)
data_Seurat <- RenameIdents(data_Seurat, new.cluster.ids)
DimPlot(data_Seurat)
data_Seurat_markers %>%rename
  group_by(cluster)%>%
  top_n(n = 100)
top10

top10 <- as.data.frame(top10)
write.csv(x = top10, row.names =TRUE, file = "top100genesbycluster.csv" )

#trajectory analysis


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


#put the tables into an assay object to plot the protein against the rnaseq
data_Seurat[["protein"]] <- CreateAssayObject(counts = data_protein)
data_Seurat <- NormalizeData(data_Seurat, assay = "protein", normalization.method = "CLR")
data_Seurat <- ScaleData(data_Seurat, assay = "protein")

# in this plot, protein (ADT) levels are on top, and RNA levels are on the bottom
#"protein_CD14","protein_CD8",                 "protein_CD45RA"
FeaturePlot(data_Seurat, features = c("protein_CD19", "protein_CD16"), min.cutoff = "q05", max.cutoff = "q95", ncol = 2)

FeaturePlot(data_Seurat, features = c("protein_CD14", "protein_CD8"), min.cutoff = "q05", max.cutoff = "q95", ncol = 2)

FeaturePlot(data_Seurat, features = "protein_CD45RA", min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
#plot the features of the protein assay on top of the rnaseq assay umap


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

p5 <-FeaturePlot(data_Seurat, features = c("protein_CD4", "protein_CD45RA",
                                           "protein_CD8", "protein_CD3", "protein_CCR7",
                                           "protein_CD19", "protein_CD16", "protein_CD56","protein_CD14"), combine =FALSE)
p5 <-lapply(X=p5, FUN = function(x)x+
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5))+
              theme(axis.title.x = element_text(size = 5))+
              theme(axis.text.y = element_text(size = 5))+
              theme(axis.text.x = element_text(size = 5)+
              theme(legend.position = "none"))
CombinePlots(plots = p5)
            

tcells <-subset(data_Seurat, idents =c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
NKcells <- subset(data_Seurat, idents = c("NK"))
bcells <- subset(data_Seurat, idents = c("B"))
monocytes <- subset(data_Seurat, idents = c("CD4+ Mono", "CD16+ Mono"))





#TRAJECTORIES
BiocManager::install(monocle3)
BiocManager::install(biocLite)
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))
BiocManager::install('SeuratWrappers')
BiocManager::install(version = "3.10")
BiocManager::install("SeuratWrappers")   
BiocManager::install("monocle")

devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")

library(Sig
library(m)
library(SingleCellExperiment)
monocle_object <- as.cell_data_set(data_seurat)
monocle_object <- cluster_cells(cds = monocle_object, reduction_method = "UMAP")
monocle_object <- learn_graph(monocle_object, use_partition = TRUE)
monocle_object <- order_cells(monocle_object,reduction_method = "UMAP")
plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = TRUE)
