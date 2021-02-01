library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "harmony_pni.rds")
pni <- readRDS(fp)
pni <- AddMetaData(pni, Idents(pni), col.name = "broad_cluster")


celltype_counts <- pni@meta.data %>%
  count(broad_cluster) %>%
  rename(total_n = "n")

cluster_order <- rev(celltype_counts$broad_cluster)

composition_plot <- pni@meta.data %>%
  count(orig.ident, broad_cluster) %>%
  left_join(celltype_counts) %>%
  mutate(prop = n / total_n) %>%
  ggplot(aes(x = factor(broad_cluster, levels=cluster_order), y = prop)) +
  geom_col(aes(fill = orig.ident)) +
  theme_bw() +
  labs(x = "", y = "Proportion of cells") +
  guides(fill=guide_legend(title="Patient")) +
  coord_flip()
# print(composition_plot)

counts_plot <- celltype_counts %>%
  ggplot(aes(factor(broad_cluster, levels=cluster_order), total_n)) +
  geom_col() +
  labs(x = "", y = "Number of cells") +
  theme_bw() +
  coord_flip()
# print(counts_plot)

combo <- composition_plot + counts_plot + plot_layout(guides = "collect")
print(combo)

cellsPerPatient_plot <- pni@meta.data %>%
  ggplot(aes(orig.ident)) +
  geom_bar(aes(fill = orig.ident)) +
  labs(x = "", y = "Number of cells") +
  guides(fill=FALSE) +
  theme_bw()
print(cellsPerPatient_plot)
