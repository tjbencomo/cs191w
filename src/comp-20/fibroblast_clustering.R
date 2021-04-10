# Description: Compare and label fibroblast subpopulations

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(tibble, include.only = "rownames_to_column")

save_results <- TRUE
data_dir <- "data"
infp <- file.path(data_dir, "seurat", "harmony-20", "harmony_combined.rds")

cells <- readRDS(infp)
cells <- subset(cells, idents = "Fibroblasts")


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

## Look for doublets
markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .75)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Likely doublet clusters
## 6 - KCs
## 7 - Myeloid
## 9 - Endothelial
notFibroClusters <- c(6, 7, 9)
notFibro <- subset(cells, idents = notFibroClusters)
fibro <- subset(cells, idents = notFibroClusters, invert = T)

## Now split by condition and analyze separately
pniFibro <- subset(fibro, condition == "pni")
jiFibro <- subset(fibro, condition == "cSCC")

## Analyze PNI fibroblasts
pniFibro <- pniFibro %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
pniFibro <- FindClusters(pniFibro, resolution = .1)

pniUmap <- DimPlot(pniFibro, reduction = "umap", label = T, label.size = 7)
print(pniUmap)
pniSampleUmap <- DimPlot(pniFibro, reduction = "umap", group.by = "orig.ident")
print(pniSampleUmap)

pni_markers <- FindAllMarkers(pniFibro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
pni_top_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Analyze non-PNI fibroblasts
jiFibro <- jiFibro %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
jiFibro <- FindClusters(jiFibro, resolution = .2)

jiUmap <- DimPlot(jiFibro, reduction = "umap", label = T, label.size = 7)
print(jiUmap)
jiSampleUmap <- DimPlot(jiFibro, reduction = "umap", group.by = "orig.ident")
print(jiSampleUmap)

ji_markers <- FindAllMarkers(jiFibro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
ji_top_markers <- ji_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

#######################################
## Finalize and save results
#######################################
notFibroMet <- notFibro@meta.data %>%
  mutate(Fibroblast_Label = case_when(
    seurat_clusters == 6 ~ "Fibroblast-KC Doublet",
    seurat_clusters == 7 ~ "Fibroblast-Myeloid Doublet",
    seurat_clusters == 9 ~ "Fibroblast-Endothelial Doublet"
  )) %>%
  mutate(keep = "No - likely doublet") %>%
  select(-starts_with("RNA"))
pniMet <- pniFibro@meta.data %>%
  mutate(Fibroblast_Label = case_when(
    seurat_clusters == 0 ~ "PNI+ Fibro C0",
    seurat_clusters == 1 ~ "Shared Fibro CXCL14+",
    seurat_clusters == 2 ~ "Shared Fibro TPPP3+",
    seurat_clusters == 3 ~ "PNI+ Fibro C3",
    seurat_clusters == 4 ~ "PNI+ Fibro C4"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))
csccMet <- jiFibro@meta.data %>%
  mutate(Fibroblast_Label = case_when(
    seurat_clusters == 0 ~ "PNI- Fibro C0",
    seurat_clusters == 1 ~ "PNI- Fibro C1",
    seurat_clusters == 2 ~ "Shared Fibro TPPP3+",
    seurat_clusters == 3 ~ "Shared Fibro CXCL14+",
    seurat_clusters == 4 ~ "PNI- Fibro C4"
  )) %>%
  mutate(keep = "Yes") %>%
  select(-starts_with("RNA"))

met <- notFibroMet %>%
  bind_rows(
    pniMet,
    csccMet
  ) %>%
  rownames_to_column(var = "cell_id")

if (save_results) {
  outdir <- file.path(data_dir, "seurat", "harmony-20")
  saveRDS(pniFibro, file.path(outdir, "pni_fibroblasts.rds"))
  saveRDS(jiFibro, file.path(outdir, "ji_fibroblasts.rds"))
  write_csv(pni_markers, file.path(outdir, "pni_fibroblast_markers.csv"))
  write_csv(ji_markers, file.path(outdir, "ji_fibroblast_markers.csv"))
  write_csv(met, file.path(outdir, "fibroblast_metadata.csv"))
}

## Exploring CAF markers from
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6972582/#ijc32193-bib-0023
plot_expression <- function(marker, pni, ji) {
  p1 <- VlnPlot(ji, marker)
  p2 <- VlnPlot(pni, marker)
  p1 + p2
}

ji_clusters <- c(
  "Ji-0",
  "Ji-1",
  "TPPP3+",
  "CXCL14+",
  "Ji-4"
)

names(ji_clusters) <- levels(jiFibro)
jiFibro <- RenameIdents(jiFibro, ji_clusters)

pni_clusters <- c(
  "Lee-0",
  "CXCL14+",
  "TPPP3+",
  "Lee-3",
  "Lee-4"
)

names(pni_clusters) <- levels(pniFibro)
pniFibro <- RenameIdents(pniFibro, pni_clusters)

## CAF Markers
plot_expression("FAP", pniFibro, jiFibro)
plot_expression("ACTA2", pniFibro, jiFibro)
plot_expression("MFAP5", pniFibro, jiFibro)
plot_expression("COL11A1", pniFibro, jiFibro)
plot_expression("TNC", pniFibro, jiFibro)
plot_expression("PDPN", pniFibro, jiFibro)
plot_expression("ITGA11", pniFibro, jiFibro)
plot_expression("CSPG4", pniFibro, jiFibro)

## Normal Fibroblast Markers
plot_expression("PDGFRA", pniFibro, jiFibro)
plot_expression("VIM", pniFibro, jiFibro)
plot_expression("S100A4", pniFibro, jiFibro)
plot_expression("POSTN", pniFibro, jiFibro)
plot_expression("COL1A1", pniFibro, jiFibro)

## Non-Fibroblast Markers
plot_expression("EPCAM", pniFibro, jiFibro)
plot_expression("CALD1", pniFibro, jiFibro)
plot_expression("SMTN", pniFibro, jiFibro)
plot_expression("PTPRC", pniFibro, jiFibro)
plot_expression("PECAM1", pniFibro, jiFibro)
