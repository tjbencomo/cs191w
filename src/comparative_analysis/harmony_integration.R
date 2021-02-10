library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(harmony)
library(future)
library(purrr, include.only = "map_dfr")
library(tibble, include.only = "rownames_to_column")

options(future.globals.maxSize = 8000 * 1024^2)
plan("multicore", workers = 2)

data_dir <- file.path("data", "ji")
data_dir <- "/scratch/users/tbencomo/cs191w/cellranger/ji"
cells <- readRDS("data/seurat/combined_raw.rds")

cells <- cells %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(vars.to.regress = "percent.mt") %>%
    RunPCA(npcs = 30)

cells <- RunHarmony(cells, "condition", max.iter.harmony = 30, max.iter.cluster = 50) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = .6)

p1 <- DimPlot(cells, reduction = "umap", label = T)
p2 <- DimPlot(cells, reduction = "umap", group.by = "condition")
p3 <- DimPlot(cells, reduction = "umap", group.by = "orig.ident")
combo_plot <- p1 + p2

ggsave("figures/harmony_umap.eps", combo_plot, width = 12, height = 6)
ggsave("figures/harmony_samples.eps", p3)

print("Finished saving images...")

get_conserved <- function(cluster){
    FindConservedMarkers(cells,
    ident.1 = cluster,
    grouping.var = "condition",
    only.pos = TRUE, logfc.threshold = .5) %>%
    rownames_to_column(var = "gene") %>%
    mutate(cluster_id = cluster)
}

cluster_ids <- 0:(length(unique(cells$seurat_clusters))-1)
markers <- map_dfr(cluster_ids, get_conserved)
top_markers <- markers %>%
  mutate(combined_avg_log2FC = (pni_avg_log2FC + cSCC_avg_log2FC) / 2) %>%
  group_by(cluster_id) %>%
  slice_max(combined_avg_log2FC, n=10) %>%
  ungroup() %>%
  select(gene, combined_avg_log2FC, cluster_id, max_pval, everything())

# markers <- FindAllMarkers(cells, only.pos = T, min.pct = .25, logfc.threshold = .5)
# top_markers <- markers %>%
#     filter(p_val_adj < .05) %>%
#     group_by(cluster) %>%
#     slice_max(avg_log2FC, n = 10)

print(cells)
print(table(cells$orig.ident))

# p1 <- DimPlot(cells, reduction = "umap", label = T)
# p2 <- DimPlot(cells, reduction = "umap", group.by = "condition")
# p3 <- DimPlot(cells, reduction = "umap", group.by = "orig.ident")
# combo_plot <- p1 + p2

# ggsave("figures/harmony_umap.eps", combo_plot, width = 12, height = 6)
# ggsave("figures/harmony_samples.eps", p3)

outdir = file.path("data", "seurat")
outfp <- file.path(outdir, "harmony_combined.rds")

saveRDS(cells, file = outfp)

write_csv(markers, file.path(outdir, "all_harmony_markers.csv"))
write_csv(top_markers, file.path(outdir, "top_harmony_markers.csv"))

