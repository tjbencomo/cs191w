# Description: Run this after performing first time harmony integration
# on combined_raw.rds. Labels broad clusters for further fine-grained clustering

library(Seurat)
library(dplyr)
library(readr)

data_dir <- "data"
infp <- file.path(data_dir, "seurat", "harmony-20", "harmony_combined.rds")
top_marker_fp <- file.path(data_dir, "seurat", "harmony-20", "top_harmony_markers.csv")
markers_fp <- file.path(data_dir, "seurat", "harmony-20", "all_harmony_markers.csv")
top_broad_fp <- file.path(data_dir, "seurat", "harmony-20", "top_broad_markers.csv")
broad_fp <- file.path(data_dir, "seurat", "harmony-20", "all_broad_markers.csv")

cells <- readRDS(infp)
top_markers <- read_csv(top_marker_fp)
markers <- read_csv(markers_fp)

p <- DimPlot(cells, reduction = "umap", label = T)
print(p)

cluster_labels <- c(
  "Myeloid cells",
  "Epithelial cells",
  "Myeloid cells",
  "Epithelial cells",
  "Myeloid cells",
  "T/NK cells",
  "T/NK cells",
  "Epithelial cells",
  "Endothelial cells",
  "Myeloid cells",
  "Epithelial cells",
  "Epithelial cells",
  "Fibroblasts",
  "Fibroblasts",
  "Melanocytes",
  "Myeloid cells",
  "Epithelial cells",
  "Myeloid cells",
  "Epithelial cells",
  "Mast cells",
  "T/NK cells",
  "B & Plasma cells",
  "Epithelial cells",
  "Epithelial cells"
)

names(cluster_labels) <- levels(cells)
cells <- RenameIdents(cells, cluster_labels)

clusterid2label <- tibble(
  cluster_id = 0:(length(cluster_labels)-1),
  label = cluster_labels
)

markers <- markers %>%
  left_join(clusterid2label)
top_markers <- top_markers %>%
  left_join(clusterid2label)

## Save Seurat object - overwrite old one
## We only updated labels
saveRDS(cells, file = infp)

write_csv(top_markers, top_broad_fp)
write_csv(markers, broad_fp)




