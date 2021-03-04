# File: figure2.R
# Description: Generate plots for figure 2 of report
# Figure 2 describes novel immune KC population

library(Seurat)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(Nebulosa)

##############################################
## Helper Functions
##############################################

make_ranks <- function(results) {
  r <- results %>%
    arrange(avg_log2FC) %>%
    pull(avg_log2FC)
  names(r) <- results %>%
    arrange(avg_log2FC) %>%
    pull(gene)
  return(r)
}

run_fgsea <- function(gene_ranking) {
  require(fgsea)
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

##############################################
## Figure Code
##############################################


data_dir <- file.path("data", "seurat", "harmony")
pni_fp <- file.path(data_dir, "pni_kcs.rds")
pni_markers_fp <- file.path(data_dir, "pni_markers.csv")
nonpni_fp <- file.path(data_dir, "nonpni_kcs.rds")
novel_sig_fp <- file.path(data_dir, "novel_signature.csv")
novel_es_fp <- file.path(data_dir, "novel_enrichment.csv")

pni <- readRDS(pni_fp)
pni <- subset(pni, idents = "KC_Unknown", invert = T)
pni_markers <- read_csv(pni_markers_fp)
novel_sig <- read_csv(novel_sig_fp)
novel_es <- read_csv(novel_es_fp)
nonpni <- readRDS(nonpni_fp)

label_order <- c("Basal", "Cycling", "Differentiating", 
                 "TSK", "IEK")
umapdf <- pni[["umap"]]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_id") %>%
  mutate(cell_label = Idents(pni)) %>%
  mutate(dataset = "PNI") %>%
  bind_rows(
    nonpni[["umap"]]@cell.embeddings %>%
      as.data.frame() %>%
      rownames_to_column(var = "cell_id") %>%
      mutate(cell_label = Idents(nonpni),
             dataset = "No PNI")
  ) %>%
  mutate(cell_label = recode(
    cell_label,
    KC_Diff = "Differentiating",
    KC_Novel = "IEK",
    KC_Cycling = "Cycling",
    KC_Basal = "Basal",
    KC_TSK = "TSK"
  )) %>%
  mutate(
    cell_label = factor(cell_label, levels = label_order)
  )

## TOP FIGURE
tsk_blue <- rgb(1, 176, 246, maxColorValue = 255)
basal_red <- rgb(248, 118, 109, maxColorValue = 255)
diff_green <- rgb(2, 191, 125, maxColorValue = 255)
immune_pink <- rgb(231, 107, 243, maxColorValue = 255)
cycling_yellow <- rgb(163, 165, 0, maxColorValue = 255)

## Draw UMAPs showing KC populations
pni_umap <- umapdf %>%
  filter(dataset == "PNI") %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cell_label)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(basal_red, cycling_yellow, diff_green, tsk_blue, immune_pink)) +
  ggtitle("PNI Keratinocytes") +
  guides(color = guide_legend(title = "")) +
  theme(text = element_text(size = 16))

nonpni_umap <- umapdf %>%
  filter(dataset == "No PNI") %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cell_label)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(basal_red, cycling_yellow, diff_green, tsk_blue)) +
  ggtitle("Non-PNI Keratinocytes") +
  guides(color = FALSE) +
  theme(text = element_text(size = 16))

topFig <- (nonpni_umap | pni_umap) + plot_layout(guides = "collect")


## MIDDLE FIGURE
## Plot composition of keratinocyte populations
## in each dataset
kc_comp_plot <- umapdf %>%
  count(dataset, cell_label) %>%
  inner_join(
    umapdf %>%
      count(dataset) %>%
      rename(n_total = n)
  ) %>%
  mutate(dataset = case_when(
    dataset == "No PNI" ~ "cSCC",
    dataset == "PNI" ~ "pSCC"
  )) %>%
  mutate(prop = n / n_total) %>%
  ggplot(aes(dataset, prop)) +
  geom_col(aes(fill = cell_label), position = 'fill') +
  scale_fill_manual(values = c(basal_red, cycling_yellow, diff_green, tsk_blue, immune_pink)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Proportion of cells") +
  guides(fill = FALSE) + 
  # guides(fill = guide_legend(title="KC Type")) +
  theme_bw() +
  theme(text = element_text(size = 16))
print(kc_comp_plot)

iek_markers <- c("IGHG1", "IGHA1", "IGLC2", "C1orf56")
x <- plot_density(pni, iek_markers, combine = FALSE)
marker_plot <- (x[[1]] + guides(color = FALSE) | x[[2]] + guides(color = FALSE)) / 
  (x[[3]] + guides(color = FALSE) | x[[4]]) + plot_layout(guides = "collect")
middleFig <- (kc_comp_plot | marker_plot) + plot_layout(width = c(1, 3))


## BOTTOM FIGURE

basal_markers <- c("KRT15", "CCL2", "COL17A1", "CXCL14", 
                   "DST", "FTH1", "MT2A", "IGFBP5", "THBS2")
cycling_markers <- c("STMN1", "HIST1H4C", "TUBA1B", "PTTG1", 
                     "HMGB2", "H2AFZ", "TOP2A", "UBE2C", "NUSAP1")
diff_markers <- c("KRTDAP", "KRT10", "KRT1", "S100A7", "SBSN", 
                  "DMKN", "KRT60", "LYPD3", "KRT6A", "CALML5")
tsk_markers <- c("MMP10", "PTHLH", "FEZ1", "IL24", "KCNMA1",
                 "INHBA", "MAGEA4", "NT5E", "LAMC2", "SLITRK6")
immune_markers <- pni_markers %>%
  filter(label == "KC_Novel") %>%
  slice_max(avg_log2FC, n = 10) %>%
  pull(gene)


## Plot heatmaps for PNI and No PNI
## No-PNI heatmap is much larger and takes longer to draw
## due to large number of cells
# heatmap_genes <- c(basal_markers, cycling_markers, diff_markers, 
#                    tsk_markers, immune_markers)
heatmap_genes <- c(diff_markers, immune_markers, cycling_markers,
                   basal_markers, tsk_markers)

## Rescale data so we can draw heatmap
pni <- ScaleData(pni, features = heatmap_genes, vars.to.regress = "percent.mt")
# nonpni <- ScaleData(nonpni, features = heatmap_genes, vars.to.regress = "percent.mt")
pni_heatmap <- DoHeatmap(pni, features = heatmap_genes, raster = TRUE)  
# nonpni_heatmap <- DoHeatmap(nonpni, features = heatmap_genes)

bottomFig <- pni_heatmap

## Save Figures Separately to combine in powerpoint
ggsave("figures/report/figure2_top.eps", topFig, width = 8.5, height = 3.6, units = "in")
ggsave("figures/report/figure2_middle.eps", middleFig, width = 8.5, height = 3.6, units = "in")
ggsave("figures/report/figure2_bottom.eps", bottomFig, width = 8.5, height = 4.5, units = "in")



# novel.markers <- FindMarkers(pni, ident.1 = "KC_Novel", logfc.threshold = .5) %>% 
#   tibble::rownames_to_column(var = "gene") %>%
#   filter(p_val_adj < .05) %>%
#   arrange(avg_log2FC )
# novel_ranks <- make_ranks(novel.markers)
# novel_es <- run_fgsea(novel_ranks)

# ## Immune Subgroup Plots
# sig_plot <- novel_sig %>%
#   mutate(label = case_when(
#     avg_log2FC > 1 ~ gene,
#     avg_log2FC < -2.5 ~ gene,
#     TRUE ~ ""
#   )) %>%
#   ggplot(aes(-log10(p_val_adj), avg_log2FC)) +
#   geom_point(alpha = .7) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   geom_label_repel(aes(label = label)) +
#   theme_bw() +
#   labs(x = "-Log10 Adjusted P-Value", y = "Log2 Fold Change") +
#   theme(text = element_text(size = 18))
# print(sig_plot)
# 
# keep_es <- c( 
#   "GO_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN", 
#   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#   "GO_EPIDERMAL_CELL_DIFFERENTIATION",
#   "GO_CORNIFICATION", 
#   "GO_PEPTIDE_CROSS_LINKING"
#   )
# 
# better_path_labels <- c(
#   GO_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN = "Circulating Immunoglobulin\nImmune Response",
#   HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon Gamma Response",
#   GO_EPIDERMAL_CELL_DIFFERENTIATION = "Epidermal Differentiation",
#   GO_CORNIFICATION = "Cornification",
#   GO_PEPTIDE_CROSS_LINKING = "Peptide Cross Linking"
#   )
# 
# novel_esplot <- novel_es %>%
#   filter(pathway %in% keep_es) %>%
#   mutate(pathway = recode(
#     pathway,
#     !!!better_path_labels
#   )) %>%
#   mutate(direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
#   mutate(direction = factor(direction, levels = c("Upregulated", "Downregulated"))) %>%
#   mutate(pathway = factor(pathway, levels = rev(better_path_labels))) %>%
#   ggplot(aes(-log10(padj), pathway)) +
#   geom_col(aes(fill = -log10(padj))) +
#   facet_wrap(~direction, ncol = 1, scales = "free") +
#   theme_bw() +
#   labs(x = "-Log10 Adjusted P-Value", y = "") +
#   theme(text = element_text(size = 18))
# print(novel_esplot)
