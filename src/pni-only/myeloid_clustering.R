# File: myeloid_clustering.R
# Description: Investigate myleoid subtypes in PNI
# data

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(stringr)

#######################################
## Helper Functions
#######################################

check_clusters <- function(cluster_markers, celltype_markers) {
  cluster_markers %>%
    filter(gene %in% celltype_markers) %>%
    select(gene, avg_log2FC, cluster) %>%
    group_by(cluster) %>%
    summarize(n_genes = n(), mean_avg_log2FC = mean(avg_log2FC)) %>%
    arrange(desc(n_genes), desc(mean_avg_log2FC))
}

#######################################
## Analysis Code
#######################################

data_dir <- "data"
infp <- file.path(data_dir, "seurat", "pni-only", "harmony_pni.rds")
save_res <- TRUE
seurat_outfp <- file.path(data_dir, "seurat", "pni-only", "myeloid.rds")
top_markers_fp <- file.path(data_dir, "seurat", "pni-only", "top_markers_myeloid.csv")
markers_fp <- file.path(data_dir, "seurat", "pni-only", "markers_myeloid.csv")

cells <- readRDS(infp)
cells <- subset(cells, idents = "Myeloid cells")

cells <- cells %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData(vars.to.regress = "percent.mt")

cells <- cells %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 20) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
cells <- FindClusters(cells, resolution = .6)

myeloid_umap <- DimPlot(cells, reduction = "umap", label=T, label.size=6)
print(myeloid_umap)

markers <- FindAllMarkers(cells, only.pos = TRUE, logfc.threshold = .25)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

macrophage_markers <- c("CD14", "CD163", "CD68")
mdsc_markers <- c("CD14", "S100A8", "S100A9", "TREM1")
cd1c_dc_markers <- c("CD1C", "CLEC10A")
clec9a_dc_markers <- c("CLEC9A", "CADM1", "XCR1")
langerhan_markers <-c("CD207", "CD1A", "S100B")
as_dc_markers <- c("AXL", "IGFBP5", "PPP1R14A")
plasmacytoid_dc_markers <- c("CLEC4C", "IL3RA")
migrating_dc_markers <- c("CCR7", "CCL19")

subtype_markers <- list(
  "Macrophage" = macrophage_markers,
  "MDSCs" = mdsc_markers,
  "CD1C DCs" = cd1c_dc_markers,
  "CLEC9A DCs" = clec9a_dc_markers,
  "Langerhan" = langerhan_markers,
  "AS DCs" = as_dc_markers,
  "Plasmacytoid DCs" = plasmacytoid_dc_markers,
  "Migrating DCs" = migrating_dc_markers
)

VlnPlot(cells, macrophage_markers) # Macrophages
VlnPlot(cells, mdsc_markers) #MSDCs
VlnPlot(cells, cd1c_dc_markers) # CD1C DCs
VlnPlot(cells, clec9a_dc_markers) # CLEC9A DCs
VlnPlot(cells, langerhan_markers) # Langerhan
VlnPlot(cells, as_dc_markers) # AS DCs
VlnPlot(cells, plasmacytoid_dc_markers) # Plasmacytoid DCs
VlnPlot(cells, migrating_dc_markers) # Migrating DCs

for(i in 1:length(subtype_markers)) {
  print(names(subtype_markers)[i])
  print(check_clusters(markers, subtype_markers[[i]]))
}

cluster_labels <- c(
  "Doublets",
  "Macrophages",
  "MDSCs",
  "Macrophages",
  "CD1C DCs",
  "CLEC9A DCs",
  "MDSCs",
  "MDSCs",
  "Migrating DCs",
  "Macrophages",
  "Plasmacytoid DCs",
  "Langerhans cells"
)

clusterid2label <- tibble(
  cluster = factor(0:(length(unique(cells$seurat_clusters))-1)),
  label = cluster_labels
)

names(cluster_labels) <- levels(cells)
cells <- RenameIdents(cells, cluster_labels)

top_markers <- top_markers %>%
  left_join(clusterid2label)
markers <- markers %>%
  left_join(clusterid2label)

if(save_res) {
  saveRDS(cells, seurat_outfp)
  write_csv(top_markers, top_markers_fp)
  write_csv(markers, markers_fp)
}
