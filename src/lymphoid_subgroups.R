library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "cca_pni.rds")
pni <- readRDS(fp)

lymphoid <- subset(pni, idents = "Lymphoid cells")

lymphoid <- lymphoid %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors() %>%
  FindClusters(resolution = 1.0)

p1 <- DimPlot(lymphoid, reduction = "umap", label = T)
print(p1)
p2 <- DimPlot(lymphoid, reduction = "umap", split.by = "orig.ident")
print(p2)

lymphoid.markers <- FindAllMarkers(lymphoid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lymphoid.markers <- lymphoid.markers %>%
  filter(p_val_adj < .05)
top_markers <- lymphoid.markers %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()


