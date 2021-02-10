# File: keratinoctye_diffexp.R
# Description: Perform differential expression analysis
# between PNI vs non-PNI tumors on each broad KC
# group as identified in Andrew's paper

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(readr)
library(muscat)
library(purrr, include.only = "map_dfr")
library(tibble, include.only = "rownames_to_column")

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "harmony", "keratinocytes.rds")
kcs <- readRDS(fp)


kcs <- kcs %>%
  RunPCA(npc = 20) %>%
  RunHarmony("condition", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony")
kcs <- FindClusters(kcs, resolution = .2)


p1 <- DimPlot(kcs, reduction = "umap", label = T)
p2 <- DimPlot(kcs, reduction = "umap", group.by = "condition")
p3 <- DimPlot(kcs, reduction = "umap", group.by = "orig.ident")
p1 + p2

get_conserved <- function(cluster){
  FindConservedMarkers(kcs,
                       ident.1 = cluster,
                       grouping.var = "condition",
                       only.pos = TRUE, logfc.threshold = .25) %>%
    rownames_to_column(var = "gene") %>%
    mutate(cluster_id = cluster)
}

cluster_ids <- 0:(length(unique(kcs$seurat_clusters))-1)
markers <- map_dfr(cluster_ids, get_conserved) %>%
  mutate(combined_avg_log2FC = (pni_avg_log2FC + cSCC_avg_log2FC) / 2)
top_markers <- markers %>%
  group_by(cluster_id) %>%
  slice_max(combined_avg_log2FC, n=10) %>%
  ungroup() %>%
  select(gene, combined_avg_log2FC, cluster_id, max_pval, everything())

basal_markers <- c("KRT15", "CCL2", "COL17A1", "CXCL14", 
                   "DST", "FTH1", "MT2A", "IGFBP5", "THBS2")
cycling_markers <- c("STMN1", "HIST1H4C", "TUBA1B", "PTTG1", 
                     "HMGB2", "H2AFZ", "TOP2A", "UBE2C", "NUSAP1")
diff_markers <- c("KRTDAP", "KRT10", "KRT1", "S100A7", "SBSN", 
                  "DMKN", "KRT60", "LYPD3", "KRT6A", "CALML5")
tsk_markers <- c("MMP10", "PTHLH", "FEZ1", "IL24", "KCNMA1",
                 "INHBA", "MAGEA4", "NT5E", "LAMC2", "SLITRK6")

markers %>% filter(gene %in% basal_markers) %>% select(gene, combined_avg_log2FC, cluster_id)
markers %>% filter(gene %in% cycling_markers) %>% select(gene, combined_avg_log2FC, cluster_id)
markers %>% filter(gene %in% diff_markers) %>% select(gene, combined_avg_log2FC, cluster_id)
markers %>% filter(gene %in% tsk_markers) %>% select(gene, combined_avg_log2FC, cluster_id)

cluster_labels <- c(
  "Differentiating",
  "Cycling",
  "Differentiating",
  "Basal",
  "Special Cycling",
  "TSK",
  "Differentiating",
  "Differentiating",
  "Differentiating"
)

names(cluster_labels) <- levels(kcs)
kcs <- RenameIdents(kcs, cluster_labels)
# rm(kcs)
# kcs.se <- as.SingleCellExperiment(kcs)

sce <- prepSCE(as.SingleCellExperiment(kcs), 
               kid = "ident", # subpopulation assignments
               gid = "condition",  # group IDs (ctrl/stim)
               sid = "orig.ident",   # sample IDs (ctrl/stim.1234)
               drop = TRUE)
sce$cluster_id <- "KC"
# sce <- sce[, sce$cluster_id %in% c("Differentiating", "Cycling")]


mm <- mmDS(sce, method = "dream",
           n_cells = 10, n_samples = 2,
           min_count = 1, min_cells = 20,
           n_threads = 3)
# mm <- mmDS(sce, method = "nbinom")
