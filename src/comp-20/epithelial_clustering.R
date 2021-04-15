# Description: Isolate epithelial cells from joint clustering. Identify epithelial
# subpopulations - keratinocytes, pilosebaceous/eccrine cells, and any misclassified
# or doublet populations. Remove non-keratinocyte cells and then separately analyze
# PNI and common SCC keratinocytes to identify keratinocyte states.


library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(tibble, include.only = "rownames_to_column")

save_results <- TRUE
data_dir <- "data"
infp <- file.path(data_dir, "seurat", "harmony-20", "harmony_combined.rds")

cells <- readRDS(infp)
cells <- subset(cells, idents = "Epithelial cells")

cells <- cells %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
cells <- FindClusters(cells, resolution = .3)

DimPlot(cells, reduction = "umap", label = T)
DimPlot(cells, reduction = "umap", group.by = "condition")
DimPlot(cells, reduction = "umap", group.by = "orig.ident")

markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .75)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## 6 - Immune doublets
## 7 - Pilosebaceous/Eccrine cells
notKCs <- subset(cells, idents = c(6, 7))
cells <- subset(cells, idents = c(6, 7), invert = T)


pni <- subset(cells, condition == "pni")
cscc <- subset(cells, condition == "cSCC")
rm(cells)

#######################################
## PNI Clustering
#######################################
pni <- pni %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
pni <- FindClusters(pni, resolution = .4)

DimPlot(pni, reduction = "umap", label = T)

pni_markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
pni_top_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Remove misclassified endothelial cells
endo <- subset(pni, idents = 10)
pni <- subset(pni, idents = 10, invert = T)

pni <- pni %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
pni <- FindClusters(pni, resolution = .4)

DimPlot(pni, reduction = "umap", label = T)

pni_markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
pni_top_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

#######################################
## cSCC Clustering
#######################################
cscc <- cscc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
cscc <- FindClusters(cscc, resolution = .4)

DimPlot(cscc, reduction = "umap", label = T)

cscc_markers <- FindAllMarkers(cscc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
cscc_top_markers <- cscc_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

#######################################
## Finalize and save results
#######################################
## Combine metadata for labeling
notKCMet <- notKCs@meta.data %>%
  mutate(Epithelial_Label = case_when(
    seurat_clusters == 7 ~ "Pilosebaceous/Eccrine cell",
    seurat_clusters == 6 ~ "Immune cell"
  )) %>%
  mutate(keep = case_when(
    seurat_clusters == 7 ~ "Yes",
    seurat_clusters == 6 ~ "No - likely doublet"
  ))
endoMet <- endo@meta.data %>%
  mutate(Epithelial_Label = "Endothelial cell",
         keep = "No - likely doublet") %>%
  select(-starts_with("RNA"))

pniMet <- pni@meta.data %>%
  mutate(Epithelial_Label = case_when(
    seurat_clusters == 0 ~ "KC_Basal",
    seurat_clusters == 1 ~ "KC_Novel-1",
    seurat_clusters == 2 ~ "KC_TSK",
    seurat_clusters %in% c(5, 7, 9) ~ "KC_Cycling",
    seurat_clusters == 4 ~ "KC_Novel-2",
    seurat_clusters %in% c(3, 6, 8) ~ "KC_Diff",
    seurat_clusters == 10 ~ "KC_Unknown-2"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))

# pniMet <- pni@meta.data %>%
#   mutate(Epithelial_Label = "Keratinocyte",
#          keep = "Yes") %>%
#   select(-starts_with("RNA"))

csccMet <- cscc@meta.data %>%
  mutate(Epithelial_Label = case_when(
    seurat_clusters %in% c(0, 3, 8, 9) ~ "KC_Diff",
    seurat_clusters %in% c(1, 6) ~ "KC_Basal",
    seurat_clusters %in% c(2, 4) ~ "KC_Cycling",
    seurat_clusters == 5 ~ "KC_TSK",
    seurat_clusters == 7 ~ "KC_Unknown-2"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))

# csccMet <- cscc@meta.data %>%
#   mutate(Epithelial_Label = "Keratinocyte",
#          keep = "Yes") %>%
#   select(-starts_with("RNA"))

met <- notKCMet %>%
  bind_rows(
    endoMet,
    pniMet,
    csccMet
  ) %>%
  rownames_to_column(var = "cell_id")

if (save_results) {
  outdir <- file.path(data_dir, "seurat", "harmony-20")
  saveRDS(pni, file.path(outdir, "pni_kcs.rds"))
  saveRDS(cscc, file.path(outdir, "ji_kcs.rds"))
  write_csv(pni_markers, file.path(outdir, "pni_kc_markers.csv"))
  write_csv(cscc_markers, file.path(outdir, "ji_kc_markers.csv"))
  write_csv(met, file.path(outdir, "epithelial_metadata.csv"))
}
