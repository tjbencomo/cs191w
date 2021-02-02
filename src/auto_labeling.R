library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(celldex)
library(SingleR)

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "cca_pni.rds")
pni <- readRDS(fp)

hpca.se <- celldex::HumanPrimaryCellAtlasData()
pred.hesc <- SingleR(test = pni@assays$RNA[, 1:16365], ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)

encode.se <- celldex::BlueprintEncodeData()
pred.encode <- SingleR(test = pni@assays$RNA[, 1:16365], ref = encode.se, assay.type.test=1,
                       labels = encode.se$label.main)

df <- data.frame(
  manual_label = Idents(pni),
  hesc_label = pred.hesc$pruned.labels,
  encode_label = pred.encode$pruned.labels
)


pni <- AddMetaData(pni, pred.hesc$pruned.labels, col.name = "hesc_label")
pni <- AddMetaData(pni, pred.encode$pruned.labels, col.name = "encode_label")
DimPlot(pni, reduction = "umap", group.by = "encode_label", label=T)
