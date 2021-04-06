# Description: Fine grained clustering of melanocytes to check for
# misclassified nerve cells

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(tibble, include.only = "rownames_to_column")


save_results <- TRUE
data_dir <- "data"
infp <- file.path(data_dir, "seurat", "harmony-20", "harmony_combined.rds")

cells <- readRDS(infp)
cells <- subset(cells, idents = "Melanocytes")

cells <- cells %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(dims =  1:30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  # RunUMAP(dims = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
cells <- FindClusters(cells, resolution = .3)

DimPlot(cells, reduction = "umap", label = T, label.size = 8)
DimPlot(cells, reduction = "umap", group.by = "condition")
DimPlot(cells, reduction = "umap", group.by = "orig.ident")

markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .75)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

c3 <- subset(cells, idents = 3)
c3 <- c3 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(dims =  1:30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  # RunUMAP(dims = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
c3 <- FindClusters(c3, resolution = .3)

DimPlot(c3, reduction = "umap", label = T, label.size = 8)
DimPlot(c3, reduction = "umap", group.by = "condition")
DimPlot(c3, reduction = "umap", group.by = "orig.ident")

c3markers <- FindAllMarkers(c3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .75)
c3top <- c3markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()
## Cluster 3 cells look like doublets - half the doublets are KCs;
## the other half is myeloid cells

met <- cells@meta.data %>%
  mutate(Melanocyte_Label = case_when(
    seurat_clusters == 3 ~ "Doublet",
    TRUE ~ "Melanocyte"
  )) %>%
  mutate(keep = case_when(
    seurat_clusters == 3 ~ "No - likely doublet",
    TRUE ~ "Yes"
  )) %>%
  select(-starts_with("RNA")) %>%
  rownames_to_column(var = "cell_id")

if (save_results) {
  outdir <- file.path(data_dir, "seurat", "harmony-20")
  write_csv(met, file.path(outdir, "melanocyte_metadata.csv"))
}
