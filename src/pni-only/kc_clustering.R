# File: kc_clustering.R
# Description: Investigate keratinocyte subtypes in PNI
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
  
cells <- readRDS(infp)
kcs <- subset(cells, idents = "Keratinocytes")
rm(cells)

kcs <- kcs %>%
  RunPCA(npc = 30) %>%
  RunHarmony("orig.ident", max.iter.harmony = 20) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
kcs <- FindClusters(kcs, resolution = .35)

kc_umap <- DimPlot(kcs, reduction = "umap", label = T, label.size = 6)
print(kc_umap)

markers <- FindAllMarkers(kcs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .25)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## Define known KC markers
basal_markers <- c("KRT15", "CCL2", "COL17A1", "CXCL14", 
                   "DST", "FTH1", "MT2A", "IGFBP5", "THBS2")
cycling_markers <- c("STMN1", "HIST1H4C", "TUBA1B", "PTTG1", 
                     "HMGB2", "H2AFZ", "TOP2A", "UBE2C", "NUSAP1")
diff_markers <- c("KRTDAP", "KRT10", "KRT1", "S100A7", "SBSN", 
                  "DMKN", "KRT60", "LYPD3", "KRT6A", "CALML5")
tsk_markers <- c("MMP10", "PTHLH", "FEZ1", "IL24", "KCNMA1",
                 "INHBA", "MAGEA4", "NT5E", "LAMC2", "SLITRK6")
celltype_markers <- list("Basal" = basal_markers, "Cycling" = cycling_markers,
                         "Differentiating" = diff_markers, "TSK" = tsk_markers)

for(i in 1:length(celltype_markers)) {
  print(names(celltype_markers)[i])
  print(check_clusters(markers, celltype_markers[[i]]))
}

## Labels
cluster_labels <- c(
  "Basal",
  "Differentiating",
  "U1",
  "Cycling",
  "Cycling",
  "Immune cells",
  "TSK",
  "Basal",
  "U2",
  "Differentiating",
  "TSK"
)
