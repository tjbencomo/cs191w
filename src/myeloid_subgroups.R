library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "cca_pni.rds")
pni <- readRDS(fp)

myeloid <- subset(pni, idents = "Myeloid cells")

myeloid <- myeloid %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors() %>%
  FindClusters(resolution = .6)

p1 <- DimPlot(myeloid, reduction = "umap", label = T)
print(p1)
p2 <- DimPlot(myeloid, reduction = "umap", split.by = "orig.ident")
print(p2)

myeloid.markers <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
myeloid.markers <- myeloid.markers %>%
  filter(p_val_adj < .05)
top_markers <- myeloid.markers %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()
