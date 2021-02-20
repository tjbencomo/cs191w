# File: myeloid_clustering.R
# Description: Identify myeloid cell types

library(Seurat)
library(harmony)
library(dplyr)
library(tibble)
library(purrr)
library(SingleR)
library(celldex)
library(readr)

save_seurat_obj <- TRUE
data_dir <- "data"
fp <- file.path(data_dir, "seurat", "harmony", "harmony_combined.rds")
cells <- readRDS(fp)

myeloid <- subset(cells, idents = "Myeloid cells")
rm(cells)

print("Aftering excluding non-myeloid cells:")
print(table(myeloid$orig.ident))

myeloid <- myeloid %>%
  RunPCA(npc = 30) %>%
  RunHarmony("condition", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
myeloid <- FindClusters(myeloid, resolution = .6)

umap_plot <- DimPlot(myeloid, reduction = "umap", label = T)
print(umap_plot)
cond_plot <- DimPlot(myeloid, reduction = "umap", group.by = "condition")
print(cond_plot)
sample_plot <- DimPlot(myeloid, reduction = "umap", group.by = "orig.ident")
print(sample_plot)


get_conserved <- function(cluster, seurat_obj){
  FindConservedMarkers(seurat_obj,
                       ident.1 = cluster,
                       grouping.var = "condition",
                       only.pos = TRUE, logfc.threshold = .25) %>%
    rownames_to_column(var = "gene") %>%
    mutate(cluster_id = cluster)
}

cluster_ids <- 0:(length(unique(myeloid$seurat_clusters))-1)
markers <- map_dfr(cluster_ids, function(x) get_conserved(x, myeloid))
markers <- markers %>%
  mutate(combined_avg_log2FC = (pni_avg_log2FC + cSCC_avg_log2FC) / 2) %>%
  select(gene, combined_avg_log2FC, cluster_id, max_pval, everything())
top_markers <- markers %>%
  mutate(combined_avg_log2FC = (pni_avg_log2FC + cSCC_avg_log2FC) / 2) %>%
  group_by(cluster_id) %>%
  slice_max(combined_avg_log2FC, n=15) %>%
  ungroup()
  
## Comaring to SingleR predictions - most of them are predicted
## B cells???? They shouldnt be in that cluster?
## Don't use these results
immune_ref <- ImmGenData()
preds <- SingleR(test = as.SingleCellExperiment(myeloid), ref = immune_ref, assay.type.test=1,
                     labels = immune_ref$label.main)
preds <- as.data.frame(preds)
print(table(preds$pruned.labels))

cluster_labels <- c(
  "Macrophages",
  "Langerhans cells",
  "CD1C DCs",
  "MSDCs",
  "Macrophages",
  "MSDCs",
  "Unknown",
  "Migrating",
  "Lymphoid like",
  "Cycling",
  "CLEC9A DCs",
  "Doublets",
  "Doublets"
)

cluster2id <- tibble(
  cluster_id = 0:(length(cluster_labels)-1),
  label = cluster_labels
)

names(cluster_labels) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, cluster_labels)

top_markers <- top_markers %>%
  left_join(cluster2id)
markers <- markers %>%
  left_join(cluster2id)

## Subdivide Migrating and Cycling
migrating <- subset(myeloid, idents = "Migrating")
cycling <- subset(myeloid, idents = "Cycling")

migrating <- migrating %>%
  RunPCA(npc = 30) %>%
  RunHarmony("condition", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
migrating <- FindClusters(migrating, resolution = .25)
migrating_umap <- DimPlot(migrating, reduction = "umap", label = T)
migrating_umap
migrating_condplot <- DimPlot(migrating, reduction = "umap", group.by = "condition")
migrating_condplot

## Migrating clusters didn't show clear subtypes - keep as migrating
migrating_markers <- map_dfr(unique(migrating$seurat_clusters), function(x) get_conserved(x, migrating))


cycling <- cycling %>%
  RunPCA(npc = 30) %>%
  RunHarmony("condition", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
cycling <- FindClusters(cycling, resolution = .25)
cycling_umap <- DimPlot(cycling, reduction = "umap", label = T)
cycling_umap
cycling_condplot <- DimPlot(cycling, reduction = "umap", group.by = "condition")
cycling_condplot

## Cycling clusters didn't show clear subtypes - keep as Cycling
cycling_markers <- map_dfr(unique(cycling$seurat_clusters), function(x) get_conserved(x, cycling))

myeloid_fp <- file.path(data_dir, "seurat", "harmony", "myeloid_cells.rds")
if (save_seurat_obj) {
  saveRDS(myeloid, myeloid_fp)
  write_csv(top_markers, "data/seurat/harmony/top_myeloid_markers.csv")
  write_csv(markers, "data/seurat/harmony/all_myeloid_markers.csv")
}
