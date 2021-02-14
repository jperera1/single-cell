library(dplyr)
library(Seurat)
library(patchwork)

#This is assuming that the given genetic dataframes are the same structure, if you wanted to 
#look at structure of a given frame, you would have to initialize outside of the function

guided_clustering <- function(datafile, scdata="abc", number_of_cells="3k", 
                              QC_RNAmin = 200, QC_RNAmax = 2500, mito_percent = 5, 
                              top_genes = 10, modularity_granularity=0.5, 
                              cluster_ids = c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                                              "NK", "DC", "Platelet")){
  #Initialization 
  project_name = paste(scdata, number_of_cells, sep ="")
  dir.create(gsub(" ", "", paste("OutputGraphs/", project_name, "/")))
  b <- gsub(" ", "", paste("GeneDatabases/", datafile, "/"))
  c <- gsub(" ", "", paste("OutputGraphs/", project_name, "/"))
  scdata.data <- Read10X(data.dir = b) #hg19 
  scdata <- CreateSeuratObject(counts = scdata.data, project = project_name, min.cells = 3, min.features = 200)
  #setwd("~/Users/jjshark2000/ProjectsR/Sampson/project_name") 
             
  #Filtering
  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")
  plot0 <- VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  saveRDS(plot0, file = gsub(" ", "", paste(c, "feature_vinplot.rds")))
  
  plot1 <- FeatureScatter(scdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(scdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  saveRDS(plot1+plot2, file = gsub(" ", "", paste(c, "feature_plot.rds")))
  scdata <- subset(scdata, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  #Normalization 
  scdata <- NormalizeData(scdata, normalization.method = "LogNormalize", scale.factor = 10000)
  scdata <- NormalizeData(scdata)
  
  #-High Feature Variability
  scdata <- FindVariableFeatures(scdata, selection.method = "vst", nfeatures = 2000)
  top <- head(VariableFeatures(scdata), top_genes)
  plot3 <- VariableFeaturePlot(scdata)
  plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, ynudge =0, xnudge=0)
  saveRDS(plot3+plot4, file = gsub(" ", "", paste(c, "high_feature_variabilityplots.rds")))

  #Scaling the Data
  all.genes <- rownames(scdata)
  scdata <- ScaleData(scdata, features = all.genes)
  scdata[["RNA"]]@scale.data #results of scaling stored here
  
  #Linear Dimension Reduction / PCA
  scdata <- RunPCA(scdata, features = VariableFeatures(object = scdata))
  print(scdata[["pca"]], dims = 1:5, nfeatures = 5) 
  VizDimLoadings(scdata, dims = 1:3, reduction = "pca") 
  plotelbow <- ElbowPlot(scdata)
  saveRDS(plotelbow, file = gsub(" ", "", paste(c, "PCA_Elbow.rds")))
  d <- DimPlot(scdata, reduction = "pca")
  saveRDS(d, file = gsub(" ", "", paste(c, "PCA_DimRedPCA.rds")))
  e<- DimHeatmap(scdata, dims = 1, cells = 500, balanced = TRUE)
  saveRDS(e, file = gsub(" ", "", paste(c, "PCA_HeatMap1.rds")))
  f<- DimHeatmap(scdata, dims = 1:15, cells = 500, balanced = TRUE)
  saveRDS(f, file = gsub(" ", "", paste(c, "PCA_HeatMap15.rds")))

  #Clustering
  h <- (readline("Please Enter Number of Principle Components: "))
  modularity_granularity=0.5
  scdata <- FindNeighbors(scdata, dims = 1:h)
  scdata <- FindClusters(scdata, resolution = modularity_granularity)
  scdata[["first_cluster"]] <- Idents(scdata)
  
  #UMAP
  scdata <- RunUMAP(scdata, dims = 1:h)
  umapplot <- DimPlot(scdata, reduction = "umap")
  saveRDS(umapplot, file = gsub(" ", "", paste(c, "umap.rds")))
  
  #Cluster/Bio-markers
  scdata.markers <- FindAllMarkers(scdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  scdata.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  #VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
  
  # you can plot raw counts as well
  #VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
  #FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
  
  library(tidyverse)
  top10 <- scdata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  DoHeatmap(scdata, features = top10$gene) + NoLegend()
  saveRDS(scdata, file = gsub(" ", "", paste(c, "topgenes_heatmap.rds")))
  
  new.cluster.ids <- cluster_ids
  names(new.cluster.ids) <- levels(pbmc)
  scdata <- RenameIdents(scdata, new.cluster.ids)
  DimPlot(scdata, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  saveRDS(scdata, file = gsub(" ", "", paste(c, "cluster_id.rds")))
}

