library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(readr)

data_dir <- "data"

sample_names <- paste("LEE0", 1:4, sep="")
sample_list <- list()

for (i in 1:length(sample_names)) {
  p <- CreateSeuratObject(
    counts = Read10X(data.dir = file.path(data_dir, sample_names[i], "filtered_feature_bc_matrix")),
    project = sample_names[i]
  )
  sample_list[[i]] <- p
}

# pni <- merge(p1, y = list(p2, p3, p4), project = "PNI")
pni <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)], add.cell.ids = sample_names)
rm(sample_list, p)
pni[["percent.mt"]] <- PercentageFeatureSet(pni, pattern = "^MT-")
pni <- subset(pni, nFeature_RNA > 200 & percent.mt < 20)

pni <- NormalizeData(pni)
pni <- FindVariableFeatures(pni, selection.method = "vst", nfeatures = 2000)
pni <- ScaleData(pni, vars.to.regress = "percent.mt")
pni <- RunPCA(pni, npcs = 30, features = VariableFeatures(object = pni))

pni <- RunHarmony(pni, "orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)

pni <- FindNeighbors(pni, reduction = "harmony", dims = 1:20)
pni <- FindClusters(pni, resolution = .1)
pni <- RunUMAP(pni, reduction = "harmony", dims = 1:20)

pni_umap <- DimPlot(pni, reduction = "umap", label = T)
print(pni_umap)
pni_sample_umap <- DimPlot(pni, reduction = "umap", split.by = "orig.ident", label = T)
print(pni_sample_umap)

markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top_markers <- markers %>%
  filter(p_val_adj < .05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

cluster_labels <- c(
  "T & NK cells",
  "Keratinocytes",
  "Endothelial cells",
  "Fibroblasts",
  "Myeloid cells",
  "Keratinocytes",
  "B & Plasma cells",
  "Melanocytes",
  "Pilosebaceous & Eccrine cells",
  "Mast cells"
)

names(cluster_labels) <- levels(pni)
pni <- RenameIdents(pni, cluster_labels)

clusterid2label <- tibble(
  cluster = factor(0:(length(unique(pni$seurat_clusters))-1)),
  label = cluster_labels
)

markers <- markers %>%
  left_join(clusterid2label)
top_markers <- top_markers %>%
  left_join(clusterid2label)


## Save Seurat object
saveRDS(pni, file = "data/seurat/pni-only/harmony_pni.rds")

write_csv(top_markers, "data/seurat/pni-only/top_markers.csv")
write_csv(markers, "data/seurat/pni-only/all_markers.csv")
# 
# ggsave("figures/harmony_umap.eps", p_tot)
# ggsave("figures/harmony_samples.eps", p_sample, width = 12, height = 6)


