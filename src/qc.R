# File: qc.R
# Description: Load PNI dataset and apply quality
# control filters before downstream analysis. QC
# includes excluding cells with high mitochondrial
# counts (assuming they are not high quality based on total depth),
# low quality cells with too few genes, and doublets.

library(Seurat)
library(BiocSingular)
library(patchwork)
# library(scater)
# library(scuttle)
# library(DropletUtils)


data_dir <- "data"

patients <- list.files(data_dir)

p1 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE01", "filtered_feature_bc_matrix")),
  project = "LEE01"
)
p2 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE02", "filtered_feature_bc_matrix")),
  project = "LEE02"
)
p3 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE03", "filtered_feature_bc_matrix")),
  project = "LEE03"
)
p4 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE04", "filtered_feature_bc_matrix")),
  project = "LEE04"
)

pni <- merge(p1, y = list(p2, p3, p4), add.cell.ids = patients, project = "PNI-SCC")
rm(p1, p2, p3, p4)
batch_info <- c(rep("batch1", sum(pni$orig.ident == "LEE01")), rep("batch2", sum(pni$orig.ident != "LEE01")))
pni <- AddMetaData(pni, metadata = batch_info, col.name = "batch")

pni[["percent.mt"]] <- PercentageFeatureSet(pni, pattern = "^MT-")


metric_plot <- VlnPlot(pni, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
depth_plot <- FeatureScatter(pni, "percent.mt", "nCount_RNA", group.by = "batch")
feature_plot <- FeatureScatter(pni, "percent.mt", "nFeature_RNA", group.by = "batch")

metric_plot / (depth_plot + feature_plot)

# ## Doublet Detection
# 
# p1.sce <- as.SingleCellExperiment(p1)
# 
# set.seed(100)
# 
# # Setting up the parameters for consistency with denoisePCA();
# # this can be changed depending on your feature selection scheme.
# dbl.dens <- computeDoubletDensity(p1.sce, subset.row=top.mam, 
#                                   d=ncol(reducedDim(sce.mam)))
# summary(dbl.dens)


# pni <- subset(pni, subset = nFeature_RNA > 200 & percent.mt < 15)

# pni <- NormalizeData(pni)
# pni <- FindVariableFeatures(pni, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(pni)
# pni <- ScaleData(pni, features = all.genes)
# pni <- RunPCA(pni, features = VariableFeatures(object = pni))
# 
# pni <- FindNeighbors(pni, dims = 1:20)
# pni <- FindClusters(pni, resolution = 0.1)
# pni <- RunUMAP(pni, dims = 1:20)
# DimPlot(pni, reduction = "umap")

# VlnPlot(big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
# s <- subset(big, subset = nFeature_RNA > 200 & percent.mt < 15)


