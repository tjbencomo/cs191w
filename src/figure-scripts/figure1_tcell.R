# File: figure1_tcell.R
# Description: Create T cell figures for Figure 1

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

data_dir <- file.path("data", "seurat", "pni-only")
infp <- file.path(data_dir, "t_cells.rds")

tcells <- readRDS(infp)
tcells <- subset(tcells, idents = "Unknown", invert = T)
tcells@meta.data$celltype <- Idents(tcells)

## T cell UMAP
tcell_umap <- DimPlot(tcells, reduction = "umap")

celltype_counts <- tcells@meta.data %>%
  count(celltype) %>%
  rename(n_cells = n) %>%
  arrange(desc(n_cells))

celltype_order <- celltype_counts %>%
  pull(celltype) %>%
  rev()


## Cell Type Composition By Patient
celltype_comp_plot <- tcells@meta.data %>%
  count(celltype, orig.ident) %>%
  left_join(celltype_counts) %>%
  mutate(celltype = factor(celltype, levels = celltype_order)) %>%
  mutate(prop = n / n_cells) %>%
  ggplot(aes(prop, celltype)) +
  geom_col(aes(fill = orig.ident)) +
  theme_bw() + 
  labs(x = "", y = "") +
  guides(fill=guide_legend(title="Patient"))


## Number of cells per celltype
celltype_count_plot <- celltype_counts %>%
  mutate(celltype = factor(celltype, levels = celltype_order)) %>%
  ggplot(aes(n_cells, celltype)) +
  geom_col() +
  theme_bw() +
  labs(x = "", y = "") +
  theme(axis.text.y = element_blank())

tcell_info_plot <- celltype_comp_plot + celltype_count_plot #+
  # plot_layout(guides = "collect")
