# File: epithelial_subgroups.R
# Description: Isolate epithelial cells and identify
# subpopulations. Based on Andrew Ji's paper, we are
# looking for 3 main cell types: 1) pilosebaceous (hair/sebum)
# cells, 2) eccrine (sweat glands) 3) keratinocytes. After
# identification the goal is to remove all non-keratinocytes
# so we can perform KC specific analyses.

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "harmony_pni.rds")
pni <- readRDS(fp)
pni <- AddMetaData(pni, Idents(pni), col.name = "broad_cluster")

epi <- subset(pni, broad_cluster == "Epithelial Cells")


epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:30)
epi <- FindClusters(epi, resolution = .6)
epi <- RunUMAP(epi, reduction = "harmony", dims = 1:30)

p1 <- DimPlot(epi, reduction = "umap", label=T)
print(p1)
p2 <- DimPlot(epi, reduction = "umap", label = T, split.by = "orig.ident")
print(p2)

epi.markers <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
epi.markers <- epi.markers %>%
  filter(p_val_adj < .05)
top_markers <- epi.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)

pilosebaceous_featplot <- FeaturePlot(epi, c("SOX9", "LHX2", "KRT15"))
print(pilosebaceous_featplot)

eccrine_featplot <- FeaturePlot(epi, c("KRT7", "KRT8", "KRT14", "KRT18", "KRT19", "MUC1", 
                                       "CEACAM5", "CEACAM6", "CEACAM1"))
print(eccrine_featplot)

epi <- AddMetaData(epi, epi$seurat_clusters, col.name = "epi_clusters")
kcs <- subset(epi, epi_clusters != 0 & epi_clusters != 8 & epi_clusters != 9 & epi_clusters != 12)

## We aren't interested in non-KC cells anymore - only save KCs
# Save Seurat object
saveRDS(kcs, file = "data/seurat/keratinocytes.rds")


ggsave("figures/epi_clusters.eps", p1)
write_csv(epi.markers, "data/markers/epithelial_subgroup_markers.csv")
