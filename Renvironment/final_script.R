library(ggplot2)
library(Seurat)
library(dplyr)
library(celda)

data = read.csv("/home/genomics/single_cell_data/Protein_umi.csv", header = TRUE,  row.names = 1)
data_rna
#head(data, n=10) 
#plug in Rata by changing the link
data_protein <- data
rm(data)
data_protein_Seurat <-data_Seurat
data_protein_seurat_markers <- data_Seurat_markers
rm(data_Seurat)
rm(data_Seurat_markers)
#this counts the number of cells that are not empty in each column (cell)
#numberofexpressed = as.data.frame(colSums(data!=0))
#we want to show he number of genes expressed in each cell, this will help us understand the quality
#ggplot(numberofexpressed) +
#  geom_histogram(
#    breaks = seq(1000,10000, by = 100),
#    aes(x=(colSums(data!=0))))
#scale_x_continuous(breaks = seq(1000,10000, by = 100))
#decontx can decontaminate data sets, it didnt work though
#mat <- as.matrix(data)
#storage.mode(mat)<-"integer"
#decontX(count = data)


#creates an object of Seurat class from count matrix from library Seurat
data_Seurat <- CreateSeuratObject(counts = data, project = "PBMC", min.cells = 3, min.features= 2)

# a lot of documentation suggests to put "^MT-" it didn't work here, "MT-" works.
#grep("MT-", rownames(data_Seurat@assays$RNA@counts), value = TRUE)
#adds the percent of mitochondiral genes in each cell, a live sample will have less than 15%          
data_Seurat[["percent.mt"]] <- PercentageFeatureSet(object = data_Seurat, pattern = "MT-")
data_Seurat[["percent.mouse"]] <- PercentageFeatureSet(object = data_Seurat, pattern = "MOUSE")

#visualise QC [metricsthe wo following checks makes sure the data is in the suerat file 
#and that the datis in the right formats 

#grep("^MT", rownames(data_Seurat[["RNA"]]), value = T)

#data_Seurat@assays$RNA@counts[1:10,1:10]

#VlnPlot(f_mouseetf1, features = c("nFeature_RNA", "nCount_RNA","percent.mouse"), ncol = 4, pt.size = 0.2)


#plot1 <- FeatureScatter(data_Seurat, feature1= "nCount_RNA", feature2="percent.mt")
#plot2 <- FeatureScatter(data_Seurat, feature1= "nCount_RNA", feature2="nFeature_RNA")
#CombinePlots(plots=list(plot1,plot2))

#We need to normalise to make it easier to visualise the expression differences
data_Seurat <-NormalizeData(data_Seurat, normalization.method = "LogNormalize")

#We want to find the most variable genes 
data_Seurat <-FindVariableFeatures(Data_seurat, selection.method = "vst", nfeatures = 2000)
#topten <-head(VariableFeatures(data_Seurat),10)
#topten
#plot1 <-VariableFeaturePlot(data_Seurat)
#plot2 <- LabelPoints(plot= plot1, points = topten, repel = TRUE)
#plot2

all.genes <-rownames(data_Seurat)
data_Seurat <- ScaleData(data_Seurat, features = all.genes)
#We want to run a PCA on the data and show a few plots 

data_protein_Seurat<-RunPCA(object = data_protein_Seurat, features = VariableFeatures(object = data_Seurat))                         
print(data_Seurat[["pca"]], dims = 1:2, nfeatures = 5)      

VizDimLoadings(data_Seurat, dims = 1:2, nfeatures = 15, reduction = "pca")

DimPlot(data_Seurat, reduction ="pca")
#DimHeatmap(data_Seurat, dims = 1, cells =500, balanced = TRUE)

#we want to determine how many clusters we have in the data
ElbowPlot(data_Seurat)

data_Seurat<-FindNeighbors(data_Seurat, dims = 1:10)
data_Seurat<-FindClusters(data_Seurat, resolution = 0.4)

head(Idents(data_Seurat), 10)

ElbowPlot(data_protein_Seurat)

data_protein_Seurat<-FindNeighbors(data_protein_Seurat, dims = 1:10)
data_protein_Seurat<-FindClusters(data_protein_Seurat, resolution = 0.2)

head(Idents(data_protein_Seurat), 10)

#UMPA is better for visualising a lot of data
data_Seurat <- RunUMAP(data_Seurat, dims = 1:10)
#when you run the below, it shouls have a graph with cells coloured by cluster

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

top10 <- human_clusters.markers %>%
  group_by(cluster)%>%
  top_n(n = 10, wt= )
top10

write.csv(top10, row.names = TRUE, file ="topeten_human_clusters.csv")
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

human_clusters <- subset(x=data_Seurat, idents = c(ident.1 = 5, ident.1 = 8, ident.1 = 10, indent.1 =11), invert = TRUE)
human_clusters <- subset(x=human_clusters, subset = percent.mouse < 5)
human_clusters <- subset(x= human_clusters, subset  = percent.mt <10)
human_clusters <- subset(x = human_clusters, subset = nFeature_RNA > 50)
human_clusters[["percent.mt"]] <- PercentageFeatureSet(object = human_clusters, pattern = "MT-")
human_clusters[["percent.mouse"]] <- PercentageFeatureSet(object = human_clusters, pattern = "MOUSE")

#plot the new graphs to check
#VlnPlot(human_clusters, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.2)
#VlnPlot(human_clusters, features = c("percent.mt","percent.mouse"), ncol = 2, pt.size = 0.2)

#renormalise the data 
human_clusters <-FindVariableFeatures(human_clusters, selection.method = "vst", nfeatures = 2000)
#topten <-head(VariableFeatures(human_clusters),10)
#topten
#plot1 <-VariableFeaturePlot(human_clusters)
#plot2 <- LabelPoints(plot= plot1, points = topten, xnudge = 0, ynudge = 0, repel =TRUE)
#plot

all.genes <-rownames(human_clusters)
human_clusters <- ScaleData(human_clusters, features = all.genes)
#rerun the PCA again to make sure the clusters are ok

human_clusters<-RunPCA(object = human_clusters, features = VariableFeatures(object = human_clusters))                         
#print(human_clusters[["pca"]], dims = 1:2, nfeatures = 5)     

#VizDimLoadings(human_clusters, dims = 1:2, nfeatures = 15, reduction = "pca")

#DimPlot(human_clusters, reduction ="pca")
#change the dims to 1,2,3 to get heatmaps for the respsective PCS
#DimHeatmap(human_clusters, dims = 3, cells =500, balanced = TRUE)

#we want to determine how many clusters we have in the data
#ElbowPlot(human_clusters)

human_clusters<-FindNeighbors(human_clusters, dims = 1:10)
human_clusters<-FindClusters(human_clusters, resolution = 0.2)

head(Idents(human_clusters), 16)

#UMPA is better for visualising a lot of dat
human_clusters <- RunUMAP(human_clusters, dims = 1:10)

#human_clusters <- subset(x= human_clusters, idents = c(ident.1 =8), invert = TRUE)

DimPlot(human_clusters, reduction = "umap")

#find markers again but just do find all this time 
human_clusters.markers <- FindAllMarkers(human_clusters, only.pos = TRUE, pin.pct = 0.25, logfc.threshold = 0.25)


x <-data_Seurat_markers %>% group_by(cluster) %>% top_n(n = 1)
FeaturePlot(data_Seurat, features = x$gene[1:4])
FeaturePlot(data_Seurat, features = x$gene[5:8])
FeaturePlot(data_Seurat, features = x$gene[9:12])

p <- FeaturePlot(human_clusters, features = c("HUMAN-CD3D",
                                              "HUMAN-CD3E",
                                              "HUMAN-CD4",
                                              "HUMAN-IL7R",
                                              "HUMAN-CD8B"), combine = FALSE)
              
#these genes are markers of t cells
                                                                                            "HUMAN-CD8A"), combine = FALSE)
p <-lapply(X=p, FUN = function(x)x+
             theme(plot.title = element_text(size = 8)) +
             theme(axis.title.y = element_text(size = 5))+
             theme(axis.title.x = element_text(size = 5))+
             theme(axis.text.y = element_text(size = 5))+
             theme(axis.text.x = element_text(size = 5))+
             theme(legend.position = "none"))
CombinePlots(plots = p)

tcells<-subset(human_clusters, features =c(ident.1 = 0, ident.1 = 5, ident.1 = 8)
                                              
p2<- FeaturePlot(human_clusters, features = c("HUMAN-CD79A", 
                                              "HUMAN-CD79B", id
                                              "HUMAN-MS4A1",
                                              "HUMAN-GNLY",
                                              "HUMAN-NKG7"), combine = FALSE)

#these genes are markers of bcells and nk cells 
p2 <-lapply(X=p2, FUN = function(x)x+
             theme(plot.title = element_text(size = 8)) +
             theme(axis.title.y = element_text(size = 5))+
             theme(axis.title.x = element_text(size = 5))+
             theme(axis.text.y = element_text(size = 5))+
             theme(axis.text.x = element_text(size = 5))+
             theme(legend.position = "none"))
CombinePlots(plots = p2)
p3 <- FeaturePlot(human_clusters, features = c("HUMAN-S100A8",
                                               "HUMAN-S100A9",
                                               "HUMAN-LYZ"), combine = FALSE)
#these genes are markers of monocytes
#p <- FeaturePlot(human_clusters, features = c("MS4A1","GNLY","CD3E","CD14", "FCER1A", "FCGR3A","LYZ","PPBP","CD8A"), combine = FALSE)
p3 <-lapply(X=p3, FUN = function(x)x+
             theme(plot.title = element_text(size = 8)) +
             theme(axis.title.y = element_text(size = 5))+
             theme(axis.title.x = element_text(size = 5))+
             theme(axis.text.y = element_text(size = 5))+
             theme(axis.text.x = element_text(size = 5))+
             theme(legend.position = "none"))
CombinePlots(plots = p3)



top10 <- human_clusters.markers %>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_logFC)
top10
#write.csv(top10, row.names = TRUE, file ="topeen_human_clusters.csv")

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

#markers:
#  Tcell_Markers = c("HUMAN-CD3E", "HUMAN-CD3D"),
#  Bcell_Markers = c("HUMAN-CD79A", "HUMAN-CD79B", "HUMAN-MS4A1"),
#  Monocyte_Markers = c("HUMAN-S100A9", "HUMAN-S100A9", "LYZ"),
#  NKcell_markers <- "HUMAN-GNLY")
#cellTypeMappings <-list(Tcells =2, Bcells =5, Monocytes =1, NKcells= 6)

tcells <- subset(human_clusters, features = c(ident.1 = 0, ident.1 = 2, ident.1 = 5, ident.1 = 8))
nkcells <- subset(human_clusters, features = c(ident.1 = 2, ident.1 =3))
bcells <- subset(human_clusters, features = c(ident.1 = 4))
monocytes <- subset(human_clusters, features = c(ident.1 = 1 ,ident.1 = 3, ident.1 = 5, ident.1 = 7))

#rename the clusters
new.cluster.ids <-c("CD4 T cells", "CD14+ monocytes","NK cells", "CD16+ monocytes", "B cells", "Monocytes", "unknown", "CD1 DCs", "MK cells", "NK cells", "pDC")

names(new.cluster.ids) <-levels(human_clusters)
human_clusters <- RenameIdents(human_clusters, new.cluster.ids)

DimPlot(human_clusters)

#plots                               
p4 <-FeaturePlot(human_clusters, features = c("HUMAN-CD3D", 
                                          "HUMAN-SELL", "HUMAN-CREM", "HUMAN-CD8A", 
                                          "HUMAN-GNLY", "HUMAN-CD79A", "HUMAN-FCGR3A", 
                                          "HUMAN-CCL2", "HUMAN-PPBP"), combine =FALSE)
p4 <-lapply(X=p4, FUN = function(x)x+
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5))+
              theme(axis.title.x = element_text(size = 5))+
              theme(axis.text.y = element_text(size = 5))+
              theme(axis.text.x = element_text(size = 5))+
              theme(legend.position = "none"))
CombinePlots(plots = p4)
#udse this one for the protein clusters

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


p5 <-FeaturePlot(data_Seurat, features = c("CD4", "CD45RO",
                                              "HUMAN-SELL", "HUMAN-CREM", "HUMAN-CD8A", 
                                              "HUMAN-GNLY", "HUMAN-CD79A", "HUMAN-FCGR3A", 
                                              "HUMAN-CCL2", "HUMAN-PPBP"), combine =FALSE)
p5 <-lapply(X=p4, FUN = function(x)x+
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5))+
              theme(axis.title.x = element_text(size = 5))+
              theme(axis.text.y = element_text(size = 5))+
              theme(axis.text.x = element_text(size = 5))+
              theme(legend.position = "none"))
CombinePlots(plots = p5)


BiocManager::install('SummarizedExperiment', force = TRUE)
library(SummarizedExperiment)
library(devtools)

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
devtools::install_github('cole-trapnell-lab/monocle3')

tcells <- subset(x= human_clusters, idents =  c("CD4 T cells", "unknown"))
NKcells <-subset(human_clusters, idents ="NK cells")
monocytes <-subset(human_clusters, idents =c("CD14+ monocytes","CD16+ monocytes", "Monocytes"))

#now we want to analyze trajectories 
#monocle only available on R4.2 

new.cluster.ids2 <-c("1", "2", "3")
names(new.cluster.ids2) <-levels(monocytes)
monocytes <- RenameIdents(monocytes, new.cluster.ids2)

expression <-t(as.matrix(srt@assays$RNA@data))
new.cluster.ids2 <-c("1", "2", "3")
names(new.cluster.ids2) <-levels(monocytes)

install.packages("SCORPIUS")  
library(SCORPIUS)
library(Seurat)
counts <-as.matrix(monocytes@assays$RNA@data)
counts<-(t(round(2^counts)))

srt <- CreateSeuratObject(counts = counts, meta.data = monocytes@meta.data)
srt <-NormalizeData(srt)

group <-as.matrix(srt@meta.data$seurat_clusters
storage.mode(group)<-"numeric"                 

space <-reduce_dimensionality(expression, dist = "spearman", ndim =3)
draw_trajectory_plot(space, progression_group = group_name, contour = TRUE)

traj <- infer_trajectory(space)
draw_trajectory_plot(space, progression_group = group_name, path = traj$path, contour = TRUE)
