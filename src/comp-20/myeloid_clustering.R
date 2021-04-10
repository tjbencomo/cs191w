# Description: Compare and label myeloid populations
## Langerhaans/CD1C DCs not observed in PNI immune cells
## Neutrophils not observed in a first pass of non-PNI immune cells
## There may be a very small subset of neutrophils in Andrew's data;
## -> would need to do more fine grained clustering in cluster-2

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(tibble, include.only = "rownames_to_column")

save_results <- TRUE
data_dir <- "data"
infp <- file.path(data_dir, "seurat", "harmony-20", "harmony_combined.rds")

cells <- readRDS(infp)
cells <- subset(cells, idents = "Myeloid cells")

cells <- cells %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
cells <- FindClusters(cells, resolution = .2)

DimPlot(cells, label = T, label.size = 6)
DimPlot(cells, group.by = "condition")
DimPlot(cells, group.by = "orig.ident")

markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .75)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## 8 = T cell doublets - see TRAC expression
## 6/7 look like B cells - label within conditions

doublets <- subset(cells, idents = 8)
cells <- subset(cells, idents = 8, invert = T)

pni <- subset(cells, condition == "pni")
ji <- subset(cells, condition == "cSCC")
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
pni <- FindClusters(pni, resolution = .2)

pniUmap <- DimPlot(pni, label = T, label.size = 6)
print(pniUmap)
pniSampleUmap <- DimPlot(pni, group.by = "orig.ident")
print(pniSampleUmap)

pni_markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
pni_top_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

pniDoublets <- subset(pni, idents = 7)
pni <- subset(pni, idents = 7, invert = T)

pni <- pni %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
pni <- FindClusters(pni, resolution = .2)

pniUmap <- DimPlot(pni, label = T, label.size = 6)
print(pniUmap)
pniSampleUmap <- DimPlot(pni, group.by = "orig.ident")
print(pniSampleUmap)

pni_markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
pni_top_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Marker Genes
## Types = Monocyte/Macrophage, B & Plasma cells, DCs, pDCs, Monocytes, Neutrophils, Cycling
## 0 - Mono/Macro
## 1 - Mono/Macro
## 2 - B & Plasma cells; specifically B cells
## 3 - B & Plasma cells; specifically plasma
## 4 - Macrophages - group with mono/macro
## 5 - DCs - specifically CLEC9A DCs
## 6 - Neutrophils
## 7 - Cycling Myeloid
## 8 - pDCs
## 9 - Monocytes - group with mono/macro


#######################################
## Ji Clustering
#######################################
ji <- ji %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
ji <- FindClusters(ji, resolution = .2)

jiUmap <- DimPlot(ji, label = T, label.size = 6)
print(jiUmap)
jiSampleUmap <- DimPlot(ji, group.by = "orig.ident")
print(jiSampleUmap)

ji_markers <- FindAllMarkers(ji, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
ji_top_markers <- ji_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Low quality doublet of KCs + high MT percent
jiDoublets <- subset(ji, idents = 9)
ji <- subset(ji, idents = 9, invert = T)

ji <- ji %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
ji <- FindClusters(ji, resolution = .2)

jiUmap <- DimPlot(ji, label = T, label.size = 6)
print(jiUmap)
jiSampleUmap <- DimPlot(ji, group.by = "orig.ident")
print(jiSampleUmap)

ji_markers <- FindAllMarkers(ji, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
ji_top_markers <- ji_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Marker Genes
## 0 - Mono/Macro
## 1 - DCs - specifically CD1C DCs; not observed in PNI
## 2 - Mono/Macro
## 3 - Monocytes
## 4 - Langerhaans; not observed in PNI
## 5 - Cycling
## 6 - pDCs
## 7 - B & Plasma cells
## 8 - DCs - specifically CLEC9A DCs

#######################################
## Finalize and save results
#######################################
doubletMet <- doublets@meta.data %>%
  mutate(Myeloid_Label = "T Cell Doublets") %>%
  mutate(keep = "No - likely doublets") %>%
  select(-starts_with("RNA")) %>%
  bind_rows(
    pniDoublets@meta.data %>%
      mutate(Myeloid_Label = "Doublets") %>%
      mutate(keep = "No - likely doublets") %>%
      select(-starts_with("RNA"))
  ) %>%
  bind_rows(
    jiDoublets@meta.data %>%
      mutate(Myeloid_Label = "Low Quality Doublets") %>%
      mutate(keep = "No - likely doublets/low quality") %>%
      select(-starts_with("RNA"))
  )

pniMet <- pni@meta.data %>%
  mutate(Myeloid_Label = case_when(
    seurat_clusters %in% c(0, 1, 4, 9) ~ "Mono/Macro",
    seurat_clusters %in% c(2, 3) ~ "B & Plasma cells",
    seurat_clusters == 5 ~ "DCs",
    seurat_clusters == 6 ~ "Neutrophils",
    seurat_clusters == 7 ~ "Cycling Myeloid",
    seurat_clusters == 8 ~ "pDCs"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))

## Marker Genes
## 0 - Mono/Macro
## 1 - DCs - specifically CD1C DCs; not observed in PNI
## 2 - Mono/Macro
## 3 - Monocytes
## 4 - Langerhaans; not observed in PNI
## 5 - Cycling
## 6 - pDCs
## 7 - B & Plasma cells
## 8 - DCs - specifically CLEC9A DCs
jiMet <- ji@meta.data %>%
  mutate(Myeloid_Label = case_when(
    seurat_clusters %in% c(0, 2, 3, 4) ~ "Mono/Macro",
    seurat_clusters %in% c(1, 8) ~ "DCs",
    seurat_clusters == 5 ~ "Cycling Myeloid",
    seurat_clusters == 6 ~ "pDCs",
    seurat_clusters == 7 ~ "B & Plasma cells",
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))

met <- doubletMet %>%
  bind_rows(pniMet, jiMet) %>%
  rownames_to_column(var = "cell_id")

if (save_results) {
  outdir <- file.path(data_dir, "seurat", "harmony-20")
  saveRDS(pni, file.path(outdir, "pni_myeloid.rds"))
  saveRDS(ji, file.path(outdir, "ji_myeloid.rds"))
  write_csv(pni_markers, file.path(outdir, "pni_myeloid_markers.csv"))
  write_csv(ji_markers, file.path(outdir, "ji_myeloid_markers.csv"))
  write_csv(met, file.path(outdir, "myeloid_metadata.csv"))
}
