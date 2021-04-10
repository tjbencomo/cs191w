# Description: Compare and label endothelial cell populations

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(tibble, include.only = "rownames_to_column")

save_results <- TRUE
data_dir <- "data"
infp <- file.path(data_dir, "seurat", "harmony-20", "harmony_combined.rds")

cells <- readRDS(infp)
cells <- subset(cells, idents = "Endothelial cells")

cells <- cells %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors()
cells <- FindClusters(cells, resolution = .3)

p1 <- DimPlot(cells, reduction = "umap", label = T, label.size = 7)
print(p1)
p2 <- DimPlot(cells, reduction = "umap", group.by = "condition")
print(p2)

markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .75)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## 6 - Immune contamination
## 8 - Keratinocytes
## 4 - On the fence about 4; keeping them due to proximity
## to known fibroblasts and fibroblast markers at top of DEG list
## Possible KC markers are also expressed in fibroblasts so hard to tell
notEndoClusters <- c(6, 8)
notEndo <- subset(cells, idents = notEndoClusters)
cells <- subset(cells, idents = notEndoClusters, invert = T)

pniEndo <- subset(cells, condition == "pni")
jiEndo <- subset(cells, condition == "cSCC")



## PNI Analysis
pniEndo <- pniEndo %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
pniEndo <- FindClusters(pniEndo, resolution = .3)

pniUmap <- DimPlot(pniEndo, reduction = "umap", label = T, label.size = 6)
print(pniUmap)
pniSampleUmap <- DimPlot(pniEndo, reduction = "umap", group.by = "orig.ident")
print(pniSampleUmap)

pni_markers <- FindAllMarkers(pniEndo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
top_pni_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Non-PNI Analysis
jiEndo <- jiEndo %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
jiEndo <- FindClusters(jiEndo, resolution = .3)

jiUmap <- DimPlot(jiEndo, reduction = "umap", label = T, label.size = 6)
print(jiUmap)
jiSampleUmap <- DimPlot(jiEndo, reduction = "umap", group.by = "orig.ident")
print(jiSampleUmap)

ji_markers <- FindAllMarkers(jiEndo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
top_ji_markers <- ji_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

#######################################
## Finalize and save results
#######################################
notEndoMet <- notEndo@meta.data %>%
  mutate(Endothelial_Label = case_when(
    seurat_clusters == 6 ~ "Myeloid Doublets",
    seurat_clusters == 8 ~ "Keratinocyte Doublets"
  )) %>%
  mutate(keep = "No - likely doublets") %>%
  select(-starts_with("RNA"))

## Shared Clusters:
## Shared-0 = Ji-0/Lee-0/Lee-3
## Shared-1 = Ji-1/Lee-6
pniMet <- pniEndo@meta.data %>%
  mutate(Endothelial_Label = case_when(
    seurat_clusters == 0 ~ "Endo Shared-0",
    seurat_clusters == 1 ~ "Endo Lee-1",
    seurat_clusters == 2 ~ "Endo Lee-2",
    seurat_clusters == 3 ~ "Endo Shared-0",
    seurat_clusters == 4 ~ "Endo Lee-4",
    seurat_clusters == 5 ~ "Endo Cycling",
    seurat_clusters == 6 ~ "Endo Shared-1"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))

jiMet <- jiEndo@meta.data %>%
  mutate(Endothelial_Label = case_when(
    seurat_clusters == 0 ~ "Endo Shared-0",
    seurat_clusters == 1 ~ "Endo Shared-1",
    seurat_clusters == 2 ~ "Endo Cycling"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))

met <- notEndoMet %>%
  bind_rows(
    pniMet,
    jiMet
  ) %>%
  rownames_to_column(var = "cell_id")

if (save_results) {
  outdir <- file.path(data_dir, "seurat", "harmony-20")
  saveRDS(pniEndo, file.path(outdir, "pni_endothelial.rds"))
  saveRDS(jiEndo, file.path(outdir, "ji_endothelial.rds"))
  write_csv(pni_markers, file.path(outdir, "pni_endothelial_markers.csv"))
  write_csv(ji_markers, file.path(outdir, "ji_endothelial_markers.csv"))
  write_csv(met, file.path(outdir, "endothelial_metadata.csv"))
}
