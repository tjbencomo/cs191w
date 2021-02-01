library(Seurat)
library(dplyr)
library(readr)

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "harmony_pni.rds")
pni <- readRDS(fp)
pni <- AddMetaData(pni, Idents(pni), col.name = "broad_cluster")

mel <- subset(pni, broad_cluster == "Melanocytes")

mel <- FindNeighbors(mel, reduction = "harmony", dims = 1:30)
mel <- FindClusters(mel, resolution = .9)
mel <- RunUMAP(mel, reduction = "harmony", dims = 1:30)

p1 <- DimPlot(mel, reduction = "umap", label=T)
print(p1)

mel.markers <- FindAllMarkers(mel, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- mel.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
