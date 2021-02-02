# File: epithelial_subgroups.R
# Description: Isolate epithelial cells and identify
# subpopulations. We are looking for keratinocytes and
# eccrine (sweat gland) or pilosebaceous (hair/sebum) cells.
# We remove non-KCs for further KC analysis at the end.

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

data_dir <- "data"
# fp <- file.path(data_dir, "seurat", "harmony_pni.rds")
fp <- file.path(data_dir, "seurat", "cca_pni.rds")
pni <- readRDS(fp)
# pni <- AddMetaData(pni, Idents(pni), col.name = "broad_cluster")

# epi <- subset(pni, broad_cluster == "Epithelial cells")
epi <- subset(pni, idents = "Epithelial cells")


epi <- RunPCA(epi)
epi <- RunUMAP(epi, reduction = "pca", dims = 1:30)
epi <- FindNeighbors(epi)
epi <- FindClusters(epi, resolution = .6)
# epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:30)
# epi <- FindClusters(epi, resolution = .6)
# epi <- RunUMAP(epi, reduction = "harmony", dims = 1:30)

p1 <- DimPlot(epi, reduction = "umap", label=T)
print(p1)
p2 <- DimPlot(epi, reduction = "umap", label = T, split.by = "orig.ident")
print(p2)

epi.markers <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
epi.markers <- epi.markers %>%
  filter(p_val_adj < .05)
top_markers <- epi.markers %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

pilosebaceous_featplot <- FeaturePlot(epi, c("SOX9", "LHX2", "KRT15"))
print(pilosebaceous_featplot)

eccrine_featplot <- FeaturePlot(epi, c("KRT7", "KRT8", "KRT14", "KRT18", "KRT19", "MUC1", 
                                       "CEACAM5", "CEACAM6", "CEACAM1"))
print(eccrine_featplot)

# epi <- AddMetaData(epi, epi$seurat_clusters, col.name = "epi_clusters")
# kcs <- subset(epi, epi_clusters != 0 & epi_clusters != 8 & epi_clusters != 9 & epi_clusters != 12)
kcs <- subset(epi, idents = 0:9)


## We aren't interested in non-KC cells anymore - only save KCs
# Save Seurat object
saveRDS(kcs, file = "data/seurat/keratinocytes.rds")


ggsave("figures/epi_clusters.eps", p1)
write_csv(epi.markers, "data/markers/all_epithelial_subgroup_markers.csv")
write_csv(top_markers, "data/markers/top_epithelial_subgroup_markers.csv")
