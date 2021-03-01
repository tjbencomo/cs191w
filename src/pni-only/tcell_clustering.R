# File: tcell_clustering.R
# Description: Investigate T/NK cell subtypes in PNI
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
seurat_outfp <- file.path(data_dir, "seurat", "pni-only", "t_cells.rds")
top_markers_fp <- file.path(data_dir, "seurat", "pni-only", "top_markers_t_cells.csv")
markers_fp <- file.path(data_dir, "seurat", "pni-only", "markers_t_cells.csv")

cells <- readRDS(infp)
tcells <- subset(cells, idents = "T & NK cells")
rm(cells)

tcells <- tcells %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 20) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
tcells <- FindClusters(tcells, resolution = 1.0)

tcell_umap <- DimPlot(tcells, reduction = "umap", label=T, label.size = 6)
print(tcell_umap)


markers <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .25)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Taken from Andrew's paper - Figure S4
## Focus on NK/Tregs/CD4/8 Exhausted, and then other cells are either
# CD4+ or CD8+
nk_markers <- c("FCER1G", "TYROBP", "XCL1", "XCL2", "AREG", 
                "KRT81", "CTSW", "KLRD1", "GNLY", "KLRC1")
treg_markers <- c("FOXP3", "IL2RA", "CD27", "TNFRSF4", "BATF", 
                  "CTLA4", "SAT1", "LAIR2")
cd4exh_markers <- c("LAG3", "PDCD1", "MAL", "CSF2", "IFNG", 
                    "GADD45G", "IFI6", "CXCL13", "IGFL2", "NMB", 
                    "CXCL13", "FKBP5", "CHN1", 
                    "TSHZ2", "BCAS3", "NR3C1", "GRAMD1A")
cd8exh_markers <- c("CCL4", "PRF1", "GZMB", "CD8A", "GZMA", 
                    "KLRC2", "KLRC1", "KRT86", "CCL3")


cell_markers <- list("NK" = nk_markers, "Treg" = treg_markers, "CD4+ Exhausted" = cd4exh_markers, 
                     "CD8+ Exhausted" = cd8exh_markers)

marker_plot <- VlnPlot(tcells, c("CD4", "CD8A", "CD8B", "NCAM1", "FCER1G", "XCL1", "FOXP3", "IL2RA"))
print(marker_plot)

## Identify cell subtypes
for(i in 1:length(cell_markers)) {
  print(names(cell_markers)[i])
  print(check_clusters(markers, cell_markers[[i]]))
}


cluster_labels <- c(
  "CD8+ Exhausted",
  "T cells",
  "CD8+",
  "Unknown",
  "Unknown",
  "Treg",
  "CD4+ Exhausted",
  "Treg",
  "NK",
  "CD8+ Exhausted"
)

clusterid2label <- tibble(
  cluster = factor(0:(length(unique(tcells$seurat_clusters))-1)),
  label = cluster_labels
)

names(cluster_labels) <- levels(tcells)
tcells <- RenameIdents(tcells, cluster_labels)

top_markers <- top_markers %>%
  left_join(clusterid2label)
markers <- markers %>%
  left_join(clusterid2label)

if(save_res) {
  saveRDS(tcells, seurat_outfp)
  write_csv(top_markers, top_markers_fp)
  write_csv(markers, markers_fp)
}
