library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readr)

data_dir <- "data"

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

pni.list <- list(p1, p2, p3, p4)
rm(p1, p2, p3, p4)

for (i in 1:length(pni.list)) {
  pni.list[[i]][["percent.mt"]] <- PercentageFeatureSet(pni.list[[i]], pattern = "^MT-")
  pni.list[[i]] <- subset(pni.list[[i]], nFeature_RNA > 200 & percent.mt < 20)
  
  pni.list[[i]] <- NormalizeData(pni.list[[i]])
  # pni.list[[i]] <- FindVariableFeatures(pni.list[[i]], selection.method = "vst", nfeatures = 2000)
  pni.list[[i]] <- FindVariableFeatures(pni.list[[i]], selection.method = "vst", nfeatures = 4000)
}

# pni.anchors <- FindIntegrationAnchors(object.list = pni.list, dims = 1:30)
# pni.integrated <- IntegrateData(anchorset = pni.anchors, dims = 1:30)
pni.anchors <- FindIntegrationAnchors(object.list = pni.list, anchor.features = 4000, dims = 1:50)
pni.integrated <- IntegrateData(anchorset = pni.anchors, dims = 1:50)


DefaultAssay(pni.integrated) <- "integrated"

pni.integrated <- ScaleData(pni.integrated)
pni.integrated <- RunPCA(pni.integrated, npcs = 50)
pni.integrated <- FindNeighbors(pni.integrated, dims = 1:50)
pni.integrated <- FindClusters(pni.integrated, resolution = 0.07)
pni.integrated <- RunUMAP(pni.integrated, reduction = "pca", dims = 1:50)
# pni.integrated <- ScaleData(pni.integrated)
# pni.integrated <- RunPCA(pni.integrated, npcs = 30)
# pni.integrated <- FindNeighbors(pni.integrated, dims = 1:30)
# pni.integrated <- FindClusters(pni.integrated, resolution = 0.1)
# pni.integrated <- RunUMAP(pni.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pni.integrated, reduction = "umap", label = T)
print(p1)
p2 <- DimPlot(pni.integrated, reduction = "umap", split.by = "orig.ident")
print(p2)
# p1 + p2

ggsave("figures/cca_umap.eps", p1)
ggsave("figures/cca_samples.eps", p2, width = 12, height = 6)

pni.markers <- FindAllMarkers(pni.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- pni.markers %>% 
  filter(p_val_adj < .05) %>%
  group_by(cluster) %>% 
  slice_max(avg_logFC, n=15) %>%
  ungroup()

cluster_labels <- c(
  "Epithelial cells",
  "Lymphoid cells",
  "Endothelial cells",
  "Fibroblasts",
  "Myeloid cells",
  "Epithelial cells",
  "Lymphoid cells",
  "Lymphoid cells",
  "Melanocytes",
  "Lymphoid cells",
  "Granulocytes"
)

names(cluster_labels) <- levels(pni.integrated)
pni.integrated <- RenameIdents(pni.integrated, cluster_labels)

clusterid2label <- tibble(
  cluster = factor(0:(length(unique(pni.integrated$seurat_clusters))-1)),
  label = cluster_labels
)

pni.markers <- pni.markers %>%
  left_join(clusterid2label)
top_markers <- top_markers %>%
  left_join(clusterid2label)

write_csv(top_markers, "data/markers/top_cca_markers.csv")
write_csv(pni.markers, "data/markers/all_cca_markers.csv")

## Save Seurat object
saveRDS(pni.integrated, file = "data/seurat/cca_pni.rds")
