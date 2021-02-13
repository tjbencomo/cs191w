# File: keratinoctye_diffexp.R
# Description: Perform differential expression analysis
# between PNI vs non-PNI tumors on each broad KC
# group as identified in Andrew's paper

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(readr)
library(muscat)
library(stringr)
library(purrr, include.only = "map_dfr")
library(tibble, include.only = "rownames_to_column")
library(limma, include.only = "makeContrasts")

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "harmony", "keratinocytes.rds")
kcs <- readRDS(fp)


kcs <- kcs %>%
  RunPCA(npc = 20) %>%
  RunHarmony("condition", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony")
kcs <- FindClusters(kcs, resolution = .2)


p1 <- DimPlot(kcs, reduction = "umap", label = T)
p2 <- DimPlot(kcs, reduction = "umap", group.by = "condition")
p3 <- DimPlot(kcs, reduction = "umap", group.by = "orig.ident")
p1 + p2

get_conserved <- function(cluster){
  FindConservedMarkers(kcs,
                       ident.1 = cluster,
                       grouping.var = "condition",
                       only.pos = TRUE, logfc.threshold = .25) %>%
    rownames_to_column(var = "gene") %>%
    mutate(cluster_id = cluster)
}

cluster_ids <- 0:(length(unique(kcs$seurat_clusters))-1)
markers <- map_dfr(cluster_ids, get_conserved) %>%
  mutate(combined_avg_log2FC = (pni_avg_log2FC + cSCC_avg_log2FC) / 2)
top_markers <- markers %>%
  group_by(cluster_id) %>%
  slice_max(combined_avg_log2FC, n=10) %>%
  ungroup() %>%
  select(gene, combined_avg_log2FC, cluster_id, max_pval, everything())

basal_markers <- c("KRT15", "CCL2", "COL17A1", "CXCL14", 
                   "DST", "FTH1", "MT2A", "IGFBP5", "THBS2")
cycling_markers <- c("STMN1", "HIST1H4C", "TUBA1B", "PTTG1", 
                     "HMGB2", "H2AFZ", "TOP2A", "UBE2C", "NUSAP1")
diff_markers <- c("KRTDAP", "KRT10", "KRT1", "S100A7", "SBSN", 
                  "DMKN", "KRT60", "LYPD3", "KRT6A", "CALML5")
tsk_markers <- c("MMP10", "PTHLH", "FEZ1", "IL24", "KCNMA1",
                 "INHBA", "MAGEA4", "NT5E", "LAMC2", "SLITRK6")

markers %>% filter(gene %in% basal_markers) %>% select(gene, combined_avg_log2FC, cluster_id)
markers %>% filter(gene %in% cycling_markers) %>% select(gene, combined_avg_log2FC, cluster_id)
markers %>% filter(gene %in% diff_markers) %>% select(gene, combined_avg_log2FC, cluster_id)
markers %>% filter(gene %in% tsk_markers) %>% select(gene, combined_avg_log2FC, cluster_id)

subtype_labels <- c(
  "0" = "Differentiating",
  "1" = "Cycling",
  "2" = "Differentiating",
  "3" = "Basal",
  "4" = "Special Cycling",
  "5" = "TSK",
  "6" = "Differentiating",
  "7" = "Differentiating",
  "8" = "Differentiating"
)

kcs$kc_subtype <- recode(kcs$seurat_clusters, !!!subtype_labels)
cluster_labels <- rep("KC", 9)

names(cluster_labels) <- levels(kcs)
kcs <- RenameIdents(kcs, cluster_labels)
# rm(kcs)
# kcs.se <- as.SingleCellExperiment(kcs)

sce <- prepSCE(as.SingleCellExperiment(kcs), 
               kid = "ident", # subpopulation assignments
               gid = "condition",  # group IDs (ctrl/stim)
               sid = "orig.ident",   # sample IDs (ctrl/stim.1234)
               drop = TRUE)
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
mm <- model.matrix(~ 0 + sce$group_id)
dimnames(mm) <- list(sce$sample_id, levels(sce$group_id))
contrast <- makeContrasts("pni-cSCC", levels = mm)




edgeR_res_all <- pbDS(pb, design = mm, contrast = contrast)$table[[1]]$KC
edgeR_res <- edgeR_res_all %>%
  filter(p_adj.loc < .05, abs(logFC) > .5) %>%
  arrange(logFC)



deseq_res_all <- pbDS(pb, design = mm, contrast = contrast, method = "DESeq2")$table[[1]]$KC
deseq_res <- deseq_res_all %>%
  filter(p_adj.loc < .05, abs(logFC) > .5) %>%
  arrange(stat)
ranks <- deseq_res$stat
names(ranks) <- deseq_res$gene

func_enrich <- function(gene_ranking) {
  reactome <- gmtPathways(file.path("data", "genesets", "c2.cp.reactome.v7.2.symbols.gmt"))
  hallmark <- gmtPathways(file.path("data", "genesets", "h.all.v7.2.symbols.gmt"))
  kegg <- gmtPathways(file.path("data", "genesets", "c2.cp.kegg.v7.2.symbols.gmt"))
  gobp <- gmtPathways(file.path("data", "genesets", "c5.go.bp.v7.2.symbols.gmt"))
  gomf <- gmtPathways(file.path("data", "genesets", "c5.go.mf.v7.2.symbols.gmt"))
  
  react_res <- fgsea(reactome, gene_ranking, eps = 0) %>%
    filter(padj < .05)
  hallmark_res <- fgsea(hallmark, gene_ranking, eps = 0) %>%
    filter(padj < .05)
  kegg_res <- fgsea(kegg, gene_ranking, eps = 0) %>%
    filter(padj < .05)
  gobp_res <- fgsea(gobp, gene_ranking, eps = 0) %>%
    filter(padj < .05)
  gomf_res <- fgsea(gomf, gene_ranking, eps = 0) %>%
    filter(padj < .05)
  
  eres <- rbind(react_res, hallmark_res, kegg_res, gobp_res, gomf_res)
}

es <- func_enrich(ranks)

s <- str_split(es$pathway, "_")
p <- sapply(s, function(x) str_c(x[2:length(x)], collapse =" "))

es_plot <- es %>%
  mutate(direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
  mutate(direction = factor(direction, levels = c("Upregulated", "Downregulated"))) %>%
  mutate(pathway = p) %>%
  group_by(direction) %>%
  slice_max(abs(NES), n = 5) %>%
  ungroup() %>%
  ggplot(aes(-log10(padj), pathway)) +
  geom_col(aes(fill = -log10(padj))) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
  facet_wrap(~direction, ncol = 1, scales = "free") +
  theme(axis.text.y = element_text(size = 10))

deseq_res_all %>%
  mutate(sig = ifelse(p_adj.loc < .05, "padj < .05", "padj > .05")) %>%
  ggplot(aes(log10(baseMean), logFC)) +
  geom_point(aes(color = sig))


scater::plotExpression(sce[, sce$cluster_id == "KC"],
               features = "EFNA1",
               x = "sample_id", colour_by = "group_id", ncol = 3) +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
df <- data.frame(cts = assay(sce)["EFNA1", ], sample_id = sce$sample_id, condition = sce$group_id) %>%
  as_tibble()
