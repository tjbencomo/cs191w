# File: pseudobulk_de.R
# Description: Perform differential expression analysis
# on  PNI vs non-PNI tumors. Perform at the broad label
# level ie KCs, myeloid cells, T cells, fibroblasts etc

library(Seurat)
library(readr)
library(dplyr)
library(muscat)
library(ggplot2)
library(limma, include.only = "makeContrasts")
library(data.table, include.only = "rbindlist")

cells_fp <- file.path("data", "seurat", "harmony", "harmony_combined.rds")
metadata_fp <- file.path("data", "combined_cell_metadata.csv")
cells <- readRDS(cells_fp)
metadata <- read_csv(metadata_fp)

print("Checking that metadata matches cell ordering")
stopifnot(all(metadata$cell_id == rownames(cells@meta.data)))
print("Check passed!")
cells@meta.data$celltype_level1 <- metadata$celltype_level1
cells@meta.data$celltype_level2 <- metadata$celltype_level2
cells@meta.data$celltype_level3 <- metadata$celltype_level3

Idents(cells) <- "celltype_level2"
## Remove these cell types as they are not of major interest, have few
## samples with enough cells to pass QC, and increase our global test power
cells2remove <- c("Doublets", "Pilosebaceous & Eccrine", "Melanocytes", 
                  "Mast cells")
cells <- subset(cells, idents = cells2remove, invert = TRUE)


sce <- prepSCE(as.SingleCellExperiment(cells), 
               kid = "celltype_level1", # subpopulation assignments
               gid = "condition",  # group IDs (ctrl/stim)
               sid = "orig.ident",   # sample IDs (ctrl/stim.1234)
               drop = TRUE)
rm(cells)

pb <- aggregateData(sce, assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

mm <- model.matrix(~ 0 + sce$group_id)
dimnames(mm) <- list(sce$sample_id, levels(sce$group_id))
contrast <- makeContrasts("pni - cSCC", levels = mm)
pbres <- pbDS(pb, design = mm, contrast = contrast)

res <- rbindlist(pbres$table[[1]])

# pbres_filtered <- lapply(pbres$table[[1]], function(u) {
#   u <- dplyr::filter(u, p_adj.loc < 0.05 & abs(logFC) > 1)
#   dplyr::arrange(u, p_adj.loc)
# })
# res <- rbindlist(pbres_filtered)


res %>%
  mutate(isHit = ifelse(p_adj.loc < .05 & abs(logFC) > 1, "hit", "ns")) %>%
  ggplot(aes(logCPM, logFC)) +
  geom_point(aes(color = isHit), alpha = .7) +
  facet_wrap(~cluster_id, scales = "free") +
  theme_bw() +
  labs(x = "log CPM Expression", y = "log Fold Change")

resfil <- res %>%
  filter(p_adj.loc < .05, abs(logFC) > 1)

