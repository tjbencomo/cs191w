library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "keratinocytes.rds")
kcs <- readRDS(fp)

kcs <- FindNeighbors(kcs, reduction = "harmony", dims = 1:30)
kcs <- FindClusters(kcs, resolution = .25)
kcs <- RunUMAP(kcs, reduction = "harmony", dims = 1:30)

p1 <- p1 <- DimPlot(kcs, reduction = "umap", label=T)
print(p1)
p2 <- DimPlot(kcs, reduction = "umap", label = T, split.by = "orig.ident")
print(p2)


kcs.markers <- FindAllMarkers(kcs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- kcs.markers %>%
  filter(p_val_adj < .05) %>%
  group_by(cluster) %>%
  slice_max(avg_logFC, n=15) %>%
  ungroup()

## Andrew's Markers
basal_genes <- c("KRT15", "CCL2", "COL17A1", "CXCL14", "DST", 
                 "CCN1", "FTH1", "MT2A", "IGFBP5", "THBS2")
basal_featplot <- FeaturePlot(kcs, basal_genes) + 
  plot_annotation(title = "Basal Marker Genes",
                  theme = theme(plot.title = element_text(hjust = 0.5)))
print(basal_featplot)

cycling_genes <- c("STMN1", "HIST1H4C", "TUBA1B", "PTTG1", "HMGB2", 
                   "H2AFZ", "TOP2A", "UBE2C", "UBE2C", "NUSAP1", "PCLAF")
cycling_featplot <- FeaturePlot(kcs, cycling_genes) +
  plot_annotation(title = "Cycling Marker Genes",
                  theme = theme(plot.title = element_text(hjust = 0.5)))
print(cycling_featplot)


diff_genes <- c("KRTDAP", "KRT10", "KRT1", "S100A7", "SBSN", "DMKN",
                "LYPD3", "KRT6A", "CALML5")
diff_featplot <- FeaturePlot(kcs, diff_genes) +
  plot_annotation(title = "Differentiating Marker Genes",
                  theme = theme(plot.title = element_text(hjust = 0.5)))
print(diff_featplot)

tsk_genes <- c("MMP10", "PTHLH", "FEZ1", "IL24", "KCNMA1", "INHBA", 
               "MAGEA4", "NT5E", "LAMC2", "SLITRK6")
tsk_featplot <- FeaturePlot(kcs, tsk_genes) +
  plot_annotation(title = "TSK Marker Genes",
                  theme = theme(plot.title = element_text(hjust = 0.5)))
print(tsk_featplot)
