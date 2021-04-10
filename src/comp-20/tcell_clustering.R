# Description: Compare and label T cell populations

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(tibble, include.only = "rownames_to_column")

save_results <- TRUE
data_dir <- "data"
infp <- file.path(data_dir, "seurat", "harmony-20", "harmony_combined.rds")

cells <- readRDS(infp)
cells <- subset(cells, idents = "T/NK cells")

cells <- cells %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
cells <- FindClusters(cells, resolution = .3)

DimPlot(cells, label = T, label.size = 6)
DimPlot(cells, group.by = "condition")
DimPlot(cells, group.by = "orig.ident")

markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .75)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## 7 - Myeloid doublets based on LYZ expression
notTCells <- subset(cells, idents = 7)
cells <- subset(cells, idents = 7, invert = T)

pni <- subset(cells, condition == "pni")
ji <- subset(cells, condition == "cSCC")


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
VlnPlot(pni, "CD4") #CD4+ T cells
VlnPlot(pni, c("CD8A", "CD8B")) # CD8+ cytotoxic t cells
VlnPlot(pni, c("CD4", "IL2RA")) #combo = T reg
VlnPlot(pni, c("KIR2DL4", "TRAC")) #lack of TRAC = NK cells

## 8 - may be KC doublet
## 2 - No subpopulation markers + lots of ribosomal genes
## Exclude both of these from further use
## Clusters:
## 0 - CD8+
## 1 - CD8+
## 2 - Low quality
## 3 - CD4+
## 4 - T reg
## 5 - NK cells
## 6 - CD8+
## 7 - CD8+ and cycling
## 8 - KC Doublets

pniLowQ <- subset(pni, idents = c(2, 8))
pni <- subset(pni, idents = c(2, 8), invert = T)

#######################################
## Non-PNI Clustering
#######################################
ji <- ji %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
ji <- FindClusters(ji, resolution = .4)

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
VlnPlot(ji, "CD4") #CD4+ T cells
VlnPlot(ji, c("CD8A", "CD8B")) # CD8+ cytotoxic t cells
VlnPlot(ji, c("CD4", "IL2RA")) #combo = T reg
VlnPlot(ji, c("KIR2DL4", "TRAC")) #lack of TRAC = NK cells

## No low quality/doublets to be found
## 0 - CD8+
## 1 - T reg
## 2 - CD8+
## 3 - NK cells
## 4 - Cycling and T reg
## 5 - CD8+ - similar to cluster-0
## 6 - CD8+ - similar to cluster-0


#######################################
## Finalize and save results
#######################################
badCellsMet <- notTCells@meta.data %>%
  mutate(TCell_Label = "Myeloid Doublets") %>%
  mutate(keep = "No - likely doublets") %>%
  select(-starts_with("RNA")) %>%
  bind_rows(
    pniLowQ@meta.data %>%
      mutate(TCell_Label = case_when(
        seurat_clusters == 2 ~ "Low Quality",
        seurat_clusters == 8 ~ "KC Doublets"
      )) %>%
      mutate(keep = "No - poor quality/doublets") %>%
      select(-starts_with("RNA"))
  )
pniMet <- pni@meta.data %>%
  mutate(TCell_Label = case_when(
    seurat_clusters %in% c(0, 1, 6, 7) ~ "CD8+",
    seurat_clusters == 3 ~ "CD4+",
    seurat_clusters == 4 ~ "T Reg",
    seurat_clusters == 7 ~ "NK"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))
jiMet <- ji@meta.data %>%
  mutate(TCell_Label = case_when(
    seurat_clusters %in% c(0, 2, 5, 6) ~ "CD8+",
    seurat_clusters %in% c(1, 4) ~ "T Reg",
    seurat_clusters == 3 ~ "NK"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))

met <- badCellsMet %>%
  bind_rows(pniMet, jiMet) %>%
  rownames_to_column(var = "cell_id")

if (save_results) {
  outdir <- file.path(data_dir, "seurat", "harmony-20")
  saveRDS(pni, file.path(outdir, "pni_tcells.rds"))
  saveRDS(ji, file.path(outdir, "ji_tcells.rds"))
  write_csv(pni_markers, file.path(outdir, "pni_tcell_markers.csv"))
  write_csv(ji_markers, file.path(outdir, "ji_tcell_markers.csv"))
  write_csv(met, file.path(outdir, "tcell_metadata.csv"))
}
