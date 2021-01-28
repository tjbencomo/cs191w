library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
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



pni <- merge(p1, y = list(p2, p3, p4), project = "PNI")
rm(p1, p2, p3, p4)
pni[["percent.mt"]] <- PercentageFeatureSet(pni, pattern = "^MT-")
pni <- subset(pni, nFeature_RNA > 200 & percent.mt < 20)

pni <- NormalizeData(pni)
pni <- FindVariableFeatures(pni, selection.method = "vst", nfeatures = 2000)
pni <- ScaleData(pni)
pni <- RunPCA(pni, npcs = 30, features = VariableFeatures(object = pni))

# no_int <- DimPlot(pni, reduction = "pca", group.by = "orig.ident")

pni <- RunHarmony(pni, "orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)

# after_int <- DimPlot(pni, reduction = "harmony", group.by = "orig.ident")

# no_int + after_int

pni <- FindNeighbors(pni, reduction = "harmony", dims = 1:30)
pni <- FindClusters(pni, resolution = .1)
pni <- RunUMAP(pni, reduction = "harmony", dims = 1:30)

p_sample <- DimPlot(pni, reduction = "umap", split.by = "orig.ident", label = T)
p_tot <- DimPlot(pni, reduction = "umap", label = T)

pni.markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cluster_labels <- c(
  "Epithelial Cells",
  "T cells/NK cells",
  "Endothelial cells",
  "Fibroblasts",
  "Myeloid cells",
  "Epithelial Cells",
  "B cells",
  "B cells",
  "Melanocytes",
  "Unknown",
  "Epithelial Cells",
  "Mast cells"
)

names(cluster_labels) <- levels(pni)
pni <- RenameIdents(pni, cluster_labels)

clusterid2label <- tibble(
  cluster = factor(0:(length(unique(pni$seurat_clusters))-1)),
  label = cluster_labels
)

pni.markers <- pni.markers %>%
  left_join(clusterid2label)

top_markers <- pni.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)


# Save Seurat object
saveRDS(pni, file = "data/seurat/harmony_pni.rds")

write_csv(top_markers, "data/markers/top_harmony_markers.csv")
write_csv(pni.markers, "data/markers/all_harmony_markers.csv")

ggsave("figures/harmony_umap.eps", p_tot)
ggsave("figures/harmony_samples.eps", p_sample, width = 12, height = 6)


