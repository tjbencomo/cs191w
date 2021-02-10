# File: label_broad_clusters.R
# Description: Label broad clusters with cell type for
# further fine grained clustering. Writes over original RDS
# object as it is the same object besides extra labeling info.

library(Seurat)
library(dplyr)
library(readr)

data_dir <- file.path("data", "seurat", "harmony")
infp <- file.path(data_dir, "harmony_combined.rds")
top_markers_fp <- file.path(data_dir, "top_harmony_markers.csv")
all_markers_fp <- file.path(data_dir, "all_harmony_markers.csv")


cells <- readRDS(infp)
top_markers <- read_csv(top_markers_fp)
all_markers <- read_csv(all_markers_fp)

p1 <- DimPlot(cells, reduction = "umap", label = T)
p1
p2 <- DimPlot(cells, reduction = "umap", group.by = "condition")
p2
p3 <- DimPlot(cells, reduction = "umap", group.by = "orig.ident")
p3


cluster_labels <- c(
  "T cells/NK cells",
  "Epithelial cells",
  "Myeloid cells",
  "Myeloid cells",
  "Myeloid cells",
  "Epithelial cells",
  "Epithelial cells",
  "Epithelial cells",
  "Myeloid cells",
  "Myeloid cells",
  "B/Plasma cells",
  "Endothelial",
  "Myeloid cells",
  "Fibroblasts",
  "Epithelial cells",
  "Myeloid cells",
  "Melanocytes",
  "Epithelial cells",
  "Fiborblasts",
  "T cells/NK cells",
  "Mast cells"
)



clusterid2label <- tibble(
  cluster_id = 0:(length(unique(cells$seurat_clusters))-1),
  label = cluster_labels
)

names(cluster_labels) <- levels(cells)
cells <- RenameIdents(cells, cluster_labels)

top_markers <- top_markers %>%
  left_join(clusterid2label)
all_markers <- all_markers %>%
  left_join(clusterid2label)

write_csv(top_markers, top_markers_fp)
write_csv(all_markers, all_markers_fp)

saveRDS(cells, file = infp)

