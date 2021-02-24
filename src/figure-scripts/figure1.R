# File: figure1.R
# Description: Create plots for figure 1 of paper

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

data_dir <- file.path("data", "seurat", "pni-only")
infp <- file.path(data_dir, "harmony_pni.rds")

cells <- readRDS(infp)
cells@meta.data$celltype <- Idents(cells)

## Broad UMAP
broad_umap <- DimPlot(cells, reduction = "umap", label = T, label.size = 6) + NoLegend()


## Data Tables
cellsPerPatient <- cells@meta.data %>%
  count(orig.ident) %>%
  rename(n_cells = n)

celltype_counts <- cells@meta.data %>%
  count(celltype) %>%
  rename(n_cells = n)
  arrange(desc(n))

celltype_order <- celltype_counts %>%
  pull(celltype) %>%
  rev()

## Number of cells per celltype
celltype_count_plot <- celltype_counts %>%
  mutate(celltype = factor(celltype, levels = celltype_order)) %>%
  ggplot(aes(n_cells, celltype)) +
  geom_col() +
  theme_bw() +
  labs(x = "", y = "") +
  theme(axis.text.y = element_blank())

## Cell Type Composition By Patient
celltype_comp_plot <- cells@meta.data %>%
  count(celltype, orig.ident) %>%
  left_join(celltype_counts) %>%
  mutate(celltype = factor(celltype, levels = celltype_order)) %>%
  mutate(prop = n / n_cells) %>%
  ggplot(aes(prop, celltype)) +
  geom_col(aes(fill = orig.ident)) +
  theme_bw() + 
  labs(x = "", y = "") +
  guides(fill=guide_legend(title="Patient"))

## Broad UMAP Infographic
broad_info_plot <- celltype_comp_plot + celltype_count_plot +
  plot_layout(guides = "collect")

## Percent Celltypes Per Patient
patient_order <- rev(c("LEE01", "LEE02", "LEE03", "LEE04"))
patient_celltype_plot <- cells@meta.data %>%
  count(orig.ident, celltype) %>%
  left_join(cellsPerPatient) %>%
  mutate(prop = n / n_cells) %>%
  mutate(orig.ident = factor(orig.ident, levels = patient_order)) %>%
  ggplot(aes(prop, orig.ident)) +
  geom_col(aes(fill = celltype)) +
  theme_bw() +
  labs(x = "", y = "")

## Cells Per Sample
patient_cellcount_plot <- cells@meta.data %>%
  ggplot(aes(orig.ident)) +
  geom_bar(aes(fill = orig.ident)) +
  theme_bw() + 
  labs(x = "", y = "Number of cells") +
  guides(fill=FALSE)

## Number of Genes Per Sample
patient_nFeature_plot <- cells@meta.data %>%
  ggplot(aes(orig.ident, nFeature_RNA)) +
  geom_violin(aes(fill = orig.ident)) +
  theme_bw() +
  labs(x = "", y = "Number of genes detected") +
  guides(fill=FALSE)

## Number of Reads Per Sample
patient_nCount_plot <- cells@meta.data %>%
  ggplot(aes(orig.ident, nCount_RNA)) +
  geom_violin(aes(fill = orig.ident)) +
  theme_bw() +
  labs(x = "", y = "Number of Reads") +
  guides(fill=FALSE)

## Percent Mitochondrial Reads Per Sample
patient_mt_plot <- cells@meta.data %>%
  ggplot(aes(orig.ident, percent.mt)) +
  geom_violin(aes(fill = orig.ident)) +
  theme_bw() +
  labs(x = "", y = "% MT Reads") +
  guides(fill=FALSE)

## QC Metrics
qc_plot <- patient_cellcount_plot + patient_nFeature_plot + 
  patient_nCount_plot + patient_mt_plot
