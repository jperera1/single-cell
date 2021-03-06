#-Initialization------------------------------
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size/sparse.size

#GetAssayData(pbmc)
#Assay 
#cell = column 
#gene = row 

#MetaData 
# y- axis cells, x-axis genotype of mouse 

#----Filtering ----------------------------------

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#pbmc[["RNA_countover2000"]] <- PercentageFeatureSet(pbmc, features = )
#how do i initialize a filter for n_countRNA etc?

#creating anew column for mitochrondrial data within the data set including "percent.mt", based on wehther it starts with MT- 

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#^not filtered yet, we will are just initailzing violin plots for those specific pieces of information

"Filter cells that have unique feature counts over 2,500 or less than 200
We filter cells that have >5% mitochondrial counts. We will filter in the following lines of code"


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.


plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#--Normalization ---------------------------------------



pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

#-High Feature Variability -----------------------------
# We choose this because prior studies have shown focusing on these tend to yield 
#better insights with these genes in downstream analysis to highlight biological signals in sc data 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, ynudge =0, xnudge=0)
plot1 + plot2

#Per gene and then look at gene distribution that are not wide spread and filtering them out 

#---Scaling the Data-----------------------------------

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc[["RNA"]]@scale.data #results of scaling stored here


#-linear dimensional reduction------------------------
#PCA - Principle Component Analysis 
#High dimension data, assumption is that the data is distributed in that multi-dimensional space
#All the points have coordinates in all the dimensions, some dimensions are responsible for more of the variance
#Finds principal components 

# Variance refers to the spread of a data set around its mean value, while a covariance refers to the measure 
#of the directional relationship between two random variables
#look at group of genes and compress into 1-d
#One subset of PC-1 is showing clear cell populations expressing one set of the genes, and the other not expressing 
# number of PC's is important for downstream dimensionality in 2-d space
# Downstream 2-d feature map is computed through the number of PC's that are selected

#each PC reflects a group of genes that were already filtered to be highly variable.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(pbmc, dims = 1:3, reduction = "pca") 
#what do these mean?
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the Dimensionality

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

"Identifying the true dimensionality of a dataset -- can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, 
exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, 
but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example, all three 
approaches yielded similar results, but we might have been justified in choosing anything between PC 7-12 as a cutoff.
We chose 10 here, but encourage users to consider the following:
- Dendritic cell and NK aficionados may recognize that genes strongly associated with PCs 12 and 13 define rare immune subsets (i.e. MZB1 is a marker for plasmacytoid DCs). However, these groups are so rare, they are difficult to distinguish from background noise for a dataset of this size without prior knowledge.
- We encourage users to repeat downstream analyses with a different number of PCs (10, 15, or even 50!). As you will observe, the results often do not differ dramatically.
- We advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs does signifcanltly and adversely affect results."

#1. Identify PCA's through Heat Maps
#2. Use the JackStraw Method 
#     - statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff
#3. Common Heuristic - Instant, To show standard deviation within top Principle Components (dimensions)


#- Clustering --------------------------- 

"construct a KNN (K-Nearest Neighbor) graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells 
based on the shared overlap in their local neighborhoods (Jaccard similarity). -> FindNeighbors()"


"To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., 
Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function.
-> FindCluster()"

"Definitions: 
- Jaccard Similarity - the intersection divided by the union of the two sets. 
- Louvain algorithm - a method to extract communities from large networks 
- SLM - 
- Modularity- one measure of the sturcture of networks or graphs, measures the strength of division of a network into modules
  - Networks with high modularity have dense connections between the nodes within modules but sparse connections between nodes 
  in different modules. Modularity is often used in optimization methods for detecting community structure in networks. However, 
  it has been shown that modularity suffers a resolution limit and, therefore, it is unable to detect small communities.
  Optimizing the standard modularity function."


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#Create different column in metadata 
pbmc[["first_cluster"]] <- Idents(pbmc)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
#These five cells belong to the clusters given 

#-Run non-linear dimensional reduction (UMAP/tSNE) -----------------------
# within 10-d space, UMAP visualizes the 10-d as a 2-d visualization 
# machine learning 

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")


#-Cluster Biomarkers (differentially expressed features)--------------------

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

library(tidyverse)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#-Assigning Cell-Type Identity------------------------------

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "../output/pbmc3k_final.rds")




