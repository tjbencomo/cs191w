library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(readr)

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

pni <- merge(p1, y = list(p2, p3, p4), project = "PNI")
rm(p1, p2, p3, p4)

pni[["percent.mt"]] <- PercentageFeatureSet(pni, pattern = "^MT-")
pni <- subset(pni, nFeature_RNA > 200 & percent.mt < 20)

pni <- NormalizeData(pni)
pni <- FindVariableFeatures(pni, selection.method = "vst", nfeatures = 2000)
pni <- ScaleData(pni)
pni <- RunPCA(pni, npcs = 30, features = VariableFeatures(object = pni))

pni <- FindNeighbors(pni, dims = 1:30)
pni <- FindClusters(pni, resolution = 0.1)
pni <- RunUMAP(pni, dims = 1:30)

p_tot <- DimPlot(pni, reduction = "umap", label = T)
p_samples <- DimPlot(pni, reduction = "umap", split.by = "orig.ident")

ggsave("figures/noint_umap.eps", p_tot)
ggsave("figures/noint_samples.eps", p_samples, width = 12, height = 6)

pni.markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- pni.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write_csv(top_markers, "data/markers/noint_markers.csv")
