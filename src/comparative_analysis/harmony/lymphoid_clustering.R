# File: lymphoid_clustering.R
# Description: Identify lymphoid cell types

library(Seurat)
library(harmony)
library(dplyr)
library(tibble)
library(purrr)
library(readr)
library(SingleR)
library(celldex)


#######################################
## Helper Functions
#######################################

get_conserved <- function(cluster, seurat_obj){
  FindConservedMarkers(seurat_obj,
                       ident.1 = cluster,
                       grouping.var = "condition",
                       only.pos = TRUE, logfc.threshold = .25) %>%
    rownames_to_column(var = "gene") %>%
    mutate(cluster_id = cluster)
}

check_clusters <- function(cluster_markers, celltype_markers) {
  cluster_markers %>%
    filter(gene %in% celltype_markers) %>%
    select(gene, combined_avg_log2FC, cluster_id) %>%
    group_by(cluster_id) %>%
    summarize(n_genes = n(), avg_log2FC = mean(combined_avg_log2FC)) %>%
    arrange(desc(n_genes), desc(avg_log2FC))
}

#######################################
## Analysis Code
#######################################

save_seurat_obj <- TRUE
data_dir <- "data"
save_fp <- file.path(data_dir, "seurat", "harmony", "T_cells.rds")
fp <- file.path(data_dir, "seurat", "harmony", "harmony_combined.rds")
cells <- readRDS(fp)

cells <- subset(cells, idents = "T cells/NK cells")

print("Aftering excluding non-lymphoid cells:")
print(table(cells$orig.ident))

cells <- cells %>%
  RunPCA(npc = 30) %>%
  RunHarmony("condition", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
cells <- FindClusters(cells, resolution = 1.5)

umap_plot <- DimPlot(cells, reduction = "umap", label = T, label.size = 6)
umap_plot
cond_plot <- DimPlot(cells, reduction = "umap", group.by = "condition")
cond_plot

markers <- map_dfr(0:(length(unique(cells$seurat_clusters))-1), function(x) get_conserved(x, cells))
markers <- markers %>%
  mutate(combined_avg_log2FC = (pni_avg_log2FC + cSCC_avg_log2FC) / 2) %>%
  select(gene, combined_avg_log2FC, cluster_id, max_pval, everything())
top_markers <- markers %>%
  mutate(combined_avg_log2FC = (pni_avg_log2FC + cSCC_avg_log2FC) / 2) %>%
  group_by(cluster_id) %>%
  slice_max(combined_avg_log2FC, n=15) %>%
  ungroup()

## Taken from Andrew's paper - Figure S4
nk_markers <- c("FCER1G", "TYROBP", "XCL1", "XCL2", "AREG", 
                "KRT81", "CTSW", "KLRD1", "GNLY", "KLRC1")
treg_markers <- c("FOXP3", "IL2RA", "CD27", "TNFRSF4", "BATF", 
                  "CTLA4", "SAT1", "LAIR2")
cd4exh_markers <- c("LAG3", "PDCD1", "MAL", "CSF2", "IFNG", 
                    "GADD45G", "IFI6", "CXCL13")
cd4preexh_markers <- c("IGFL2", "NMB", "CXCL13", "FKBP5", "CHN1", 
                       "TSHZ2", "BCAS3", "NR3C1", "GRAMD1A")
cd8exh_markers <- c("CCL4", "PRF1", "GZMB", "CD8A", "GZMA", 
                    "KLRC2", "KLRC1", "KRT86", "CCL3")
cd8temra_markers <- c("FGFBP2", "GZMH", "KLRG1", "NKG7", "PLAC8", 
                      "FCGR3A", "PLEK", "CCL4")
cd8tem_markers <- c("GZMK", "CCL5", "CD8B", "CMC1", "SAMD3", 
                    "GIMAP7", "GIMAP4")
cd4rgcc_markers <- c("IL7R", "ANXA1", "NBAS", "CLU", "LMNA", "RGCC", 
                     "PLIN2", "ZNF331", "CREM", "SLC2A3", "KLRB1", 
                     "ZFP36L2")
cd8naive_markers <- c("ANXA1", "NBAS", "CLU", "BCLAS3", "CD8A", "CD8B")
cd4naive_markers <- c("CCR7", "SELL", "IL7R", "ANXA1", "CD4")

cell_markers <- list("NK" = nk_markers, "Treg" = treg_markers, "CD4+ Exhausted" = cd4exh_markers, 
                     "CD4+ Pre-Exhausted" = cd4preexh_markers,"CD8+ Exhausted" = cd8exh_markers, 
                     "CD8+ T_EMRA" = cd8temra_markers, "CD8+ T_EM" = cd8tem_markers, 
                     "CD4+ RGCC" = cd4rgcc_markers, "CD8+ Naive" = cd8naive_markers, 
                     "CD4+ Naive" = cd4naive_markers)

## Identify cell subtypes
for(i in 1:length(cell_markers)) {
  print(names(cell_markers)[i])
  print(check_clusters(markers, cell_markers[[i]]))
}
## Note 11/16 exhibit KRT14/16 that is specific to epithelial cells
## These are likely doublets and should be removed

cluster_labels <- c(
  "CD8+ T_EM",
  "CD4+ Exhausted",
  "CD8+ Exhausted",
  "CD8+ Naive",
  "CD4+ Pre-Exhausted",
  "NK",
  "Treg",
  "Treg",
  "CD4+ Naive",
  "Cycling",
  "CD4+ Naive",
  "Doublets",
  "CD8+ T_EMRA",
  "CD8+ Exhausted",
  "CD8+ Naive",
  "NK",
  "Doublets"
)

cluster2id <- tibble(
  cluster_id = 0:(length(cluster_labels)-1),
  label = cluster_labels
)

names(cluster_labels) <- levels(cells)
cells <- RenameIdents(cells, cluster_labels)


top_markers <- top_markers %>%
  left_join(cluster2id)
markers <- markers %>%
  left_join(cluster2id)

if (save_seurat_obj) {
  saveRDS(cells, save_fp)
  write_csv(top_markers, "data/seurat/harmony/top_Tcell_markers.csv")
  write_csv(markers, "data/seurat/harmony/all_Tcell_markers.csv")
}
