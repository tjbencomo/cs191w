library(Seurat)
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

epi.markers <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- epi.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)

write_csv(epi.markers, "data/markers/epithelial_subgroup_markers.csv")
