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

make_ranks <- function(results, subgroup) {
  de_ranks <- results %>%
    filter(cluster_id == subgroup) %>%
    arrange(logFC) %>%
    pull(logFC)
  names(de_ranks) <- results %>%
    filter(cluster_id == subgroup) %>%
    arrange(logFC) %>%
    pull(gene)
  return(de_ranks)
}


func_enrich <- function(gene_ranking, alpha_thres = .05) {
  require(fgsea)
  reactome <- gmtPathways(file.path("data", "genesets", "c2.cp.reactome.v7.2.symbols.gmt"))
  hallmark <- gmtPathways(file.path("data", "genesets", "h.all.v7.2.symbols.gmt"))
  kegg <- gmtPathways(file.path("data", "genesets", "c2.cp.kegg.v7.2.symbols.gmt"))
  gobp <- gmtPathways(file.path("data", "genesets", "c5.go.bp.v7.2.symbols.gmt"))
  gomf <- gmtPathways(file.path("data", "genesets", "c5.go.mf.v7.2.symbols.gmt"))
  
  react_res <- fgsea(reactome, gene_ranking, eps = 0) %>%
    filter(padj < alpha_thres)
  hallmark_res <- fgsea(hallmark, gene_ranking, eps = 0) %>%
    filter(padj < alpha_thres)
  kegg_res <- fgsea(kegg, gene_ranking, eps = 0) %>%
    filter(padj < alpha_thres)
  gobp_res <- fgsea(gobp, gene_ranking, eps = 0) %>%
    filter(padj < alpha_thres)
  gomf_res <- fgsea(gomf, gene_ranking, eps = 0) %>%
    filter(padj < alpha_thres)
  
  eres <- bind_rows(react_res, hallmark_res, kegg_res, gobp_res, gomf_res)
}

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
               )
rm(cells)

pb <- aggregateData(sce, assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

mm <- model.matrix(~ 0 + sce$group_id)
dimnames(mm) <- list(sce$sample_id, c(levels(sce$group_id)))
contrast <- makeContrasts("pni - cSCC", levels = mm)
pbres <- pbDS(pb, design = mm, contrast = contrast)


res <- rbindlist(pbres$table[[1]])


ma_plot <- res %>%
  mutate(isHit = ifelse(p_adj.loc < .01 & abs(logFC) > 1, "hit", "ns")) %>%
  ggplot(aes(logCPM, logFC)) +
  geom_point(aes(color = isHit), alpha = .7) +
  facet_wrap(~cluster_id, scales = "free_x") +
  theme_bw() +
  labs(x = "log CPM Expression", y = "log Fold Change")
print(ma_plot)

resfil <- res %>%
  filter(p_adj.loc < .01, abs(logFC) > 1)

kc_ranks <- make_ranks(resfil, "Epithelial cells")
kc_es <- func_enrich(kc_ranks)

myeloid_ranks <- make_ranks(resfil, "Myeloid cells")
myeloid_es <- func_enrich(myeloid_ranks)

tcell_ranks <- make_ranks(resfil, "T cells/NK cells")
tcell_es <- func_enrich(tcell_ranks)

outfp <- file.path("data/seurat/harmony/broad_de_results.csv")
write_csv(resfil, outfp)
