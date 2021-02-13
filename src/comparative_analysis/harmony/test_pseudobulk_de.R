library(muscat)
library(limma, include.only = "makeContrasts")
library(DESeq2)

fp <- "data/seurat/harmony/kcs_sce.rds"
sce <- readRDS(fp)

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
mm <- model.matrix(~ 0 + sce$group_id)
dimnames(mm) <- list(sce$sample_id, levels(sce$group_id))
contrast <- makeContrasts("pni-cSCC", levels = mm)

# res <- pbDS(pb, design = mm, contrast = contrast, method = "DESeq2")
res <- pbDS(pb, design = ~group_id, method = "DESeq2")
tbl <- res$table[[1]]
muscat_res <- tbl$KC
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

# scater::plotExpression(sce[, sce$cluster_id == "KC"],
#                features = tbl_fil$`KC`$gene[13585],
#                x = "sample_id", colour_by = "group_id", ncol = 3) +
#   guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# kcs <- readRDS("data/seurat/harmony/keratinocytes.rds")


dds <- DESeqDataSet(pb, design = ~ group_id)
dds <- DESeq(dds)

mle <- results(dds, alpha = .05) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene")
map_res <- lfcShrink(dds, coef = "group_id_pni_vs_cSCC") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene")
