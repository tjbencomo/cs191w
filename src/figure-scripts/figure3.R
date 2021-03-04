# File: figure3.R
# Description: Generate plots for figure 3 of report

library(Seurat)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(fgsea)
library(patchwork)
library(ggrepel)
library(stringr)


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

plotES <- function(es_res, db_fil = NULL, n_terms = 5) {
  s <- str_split(es_res$pathway, "_")
  db <- sapply(s, function(x) x[1])
  p <- sapply(s, function(x) str_c(x[2:length(x)], collapse =" "))
  es_res$db <- db
  es_res <- es_res %>%
    mutate(direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
    mutate(direction = factor(direction, levels = c("Upregulated", "Downregulated"))) %>%
    mutate(pathway = p)
  
  if (!is.null(db_fil)) {
    if (length(db_fil) > 1) {
      es_res <- es_res %>%
        filter(db %in% db_fil)
    } else {
      es_res <- es_res %>%
        filter(db == db_fil)
    }
  }
  
  es_res %>%
    group_by(direction) %>%
    slice_max(abs(NES), n = n_terms) %>%
    ungroup() %>%
    ggplot(aes(-log10(padj), reorder(pathway, -log10(padj)))) +
    geom_col(aes(fill = -log10(padj))) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
    facet_wrap(~direction, ncol = 1, scales = "free") +
    theme(axis.text.y = element_text(size = 10)) +
    labs(x = "-log10(Adjusted P-Value)", y = "")
}


data_dir <- file.path("data", "seurat", "harmony")
infp <- file.path(data_dir, "keratinocytes.rds")
basal_de_fp <- file.path(data_dir, "basal_de_res.csv")
diff_de_fp <- file.path(data_dir, "diff_de_res.csv")
cycling_de_fp <- file.path(data_dir, "cycling_de_res.csv")
tsk_de_fp <- file.path(data_dir, "tsk_de_res.csv")

kcs <- readRDS(infp)
kcs$celltype <- Idents(kcs)
kcs <- subset(kcs, idents = "KC_Unknown", invert = T)
myc_path <- gmtPathways("data/genesets/HALLMARK_MYC_TARGETS.gmt")
kegg_p53_path <- gmtPathways("data/genesets/KEGG_P53_SIGNALING_PATHWAY.gmt")
hall_p53_path <- gmtPathways("data/genesets/HALLMARK_P53_PATHWAY.gmt")
emt_path <- gmtPathways("data/genesets/HALLMARK_EMT.gmt")
paths <- c(myc_path, kegg_p53_path, hall_p53_path, emt_path)

basal_de <- read_csv(basal_de_fp)
cycling_de <- read_csv(cycling_de_fp)
diff_de <- read_csv(diff_de_fp)
tsk_de <- read_csv(tsk_de_fp)

kcs <- AddModuleScore(kcs, features = paths, name = c("MYC_TARGETS", "KEGG_P53_PATHWAY", "HALLMARK_P53_PATHWAY", "EMT_SCORE"))

kc_order <- c("Basal", "Differentiating", "TSK")
myc_plot <- kcs@meta.data %>%
  filter(!celltype %in% c("KC_Cycling", "KC_Novel")) %>%
  mutate(celltype = case_when(
    celltype == "KC_Basal" ~ "Basal",
    celltype == "KC_Diff" ~ "Differentiating",
    celltype == "KC_TSK" ~ "TSK"
  )) %>%
  mutate(celltype = factor(celltype, levels = kc_order)) %>%
  ggplot(aes(MYC_TARGETS1)) +
  geom_density(aes(fill = condition), alpha=.7) +
  theme_bw() +
  labs(x = "Myc Target Activity", y = "Density") +
  scale_fill_discrete(name = "", labels=c("cSCC", "pSCC")) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_text(size = 6))
print(myc_plot)

p53_plot <- kcs@meta.data %>%
  filter(!celltype %in% c("KC_Cycling", "KC_Novel")) %>%
  mutate(celltype = case_when(
    celltype == "KC_Basal" ~ "Basal",
    celltype == "KC_Diff" ~ "Differentiating",
    celltype == "KC_TSK" ~ "TSK"
  )) %>%
  mutate(celltype = factor(celltype, levels = kc_order)) %>%
  ggplot(aes(KEGG_P53_PATHWAY2)) +
  geom_density(aes(fill = condition), alpha=.7) +
  theme_bw() +
  labs(x = "p53 Pathway Activity", y = "Density") +
  scale_fill_discrete(name = "", labels=c("cSCC", "pSCC")) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_text(size = 6))
print(p53_plot)

volc_plot <- tsk_de %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(sig = ifelse(p_val_adj < .05, "padj < .05", "ns")) %>%
  mutate(lbl = ifelse(abs(avg_log2FC) > 2.6 & p_val_adj < .05, gene, "")) %>%
  ggplot(aes(avg_log2FC, logp, label = lbl)) +
  geom_point(aes(color = sig)) +
  geom_label_repel(size = 2, min.segment.length = .01) +
  theme_bw() +
  labs(x = "Log2 Fold Change", y = "-log10(Adjusted P-Value)") +
  guides(color = guide_legend(title = ""))
print(volc_plot)

tsk_ranks <- make_ranks(tsk_de %>% filter(p_val_adj < .05))
tsk_es <- run_fgsea(tsk_ranks)

keep_paths <- c("HALLMARK_P53_PATHWAY", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE", "GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
                "GO_ADAPTIVE_IMMUNE_RESPONSE", "GO_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY",
                "GO_KERATINIZATION", "GO_CORNIFICATION", "GO_KERATINOCYTE_DIFFERENTIATION")
gsea_plot <- tsk_es %>%
  filter(pathway %in% keep_paths) %>%
  plotES() +
  theme_bw() +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(size = 6))


# totalFig <- (volc_plot / p53_plot / myc_plot | gsea_plot) + plot_layout(guides = "collect")
totalFig <- (gsea_plot | p53_plot / myc_plot) + plot_layout(guides = "collect") & theme(legend.position = 'top')

ggsave("figures/report/figure3.png", totalFig, width = 9.5, height = 5.5, unit = "in")

# myc_plot <- kcs@meta.data %>%
#   filter(celltype == "KC_Diff") %>%
#   mutate(celltype = "Differentiating KCs") %>%
#   ggplot(aes(MYC_TARGETS1)) +
#   geom_density(aes(fill = condition), alpha=.7) +
#   theme_bw() +
#   labs(x = "MYC TARGETS ACTIVITY", y = "Density") +
#   scale_fill_discrete(name = "", labels=c("Non-PNI", "PNI")) +
#   facet_wrap(~celltype)
# 
# p53_plot <- kcs@meta.data %>%
#   filter(celltype == "KC_Diff") %>%
#   ggplot(aes(P53_PATHWAY2)) +
#   geom_density(aes(fill = condition), alpha=.7) +
#   theme_bw() +
#   labs(x = "P53 PATHWAY ACTIVITY", y = "Density") +
#   scale_fill_discrete(name = "", labels=c("Non-PNI", "PNI"))
# 
# activity_plots <- myc_plot / p53_plot + plot_layout(guides = "collect")
# print(activity_plots)

## DE Analysis of KC Subtypes
# de <- basal_de %>%
#   bind_rows(cycling_de) %>%
#   bind_rows(diff_de) %>%
#   bind_rows(tsk_de) %>%
#   mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
#   mutate(id = str_c(gene, direction, sep = "_"))
# 
# de %>%
#   mutate(sig = ifelse(p_val_adj < .05, "padj < .05", "ns")) %>%
#   mutate(sig = factor(sig, levels = c("padj < .05", "ns"))) %>%
#   mutate(lbl = case_when(
#     subgroup == "Basal" & (avg_log2FC > 2 | avg_log2FC < -2) & p_val_adj < .05 ~ gene,
#     subgroup == "Cycling" & (avg_log2FC > 2 | avg_log2FC < -2) & p_val_adj < .05 ~ gene,
#     subgroup == "Differentiating" & (avg_log2FC > 4 | avg_log2FC < -2) & p_val_adj < .05 ~ gene,
#     subgroup == "TSK" & (avg_log2FC > 2.5 | avg_log2FC < -2.3) & p_val_adj < .05 ~ gene,
#     TRUE ~ ""
#   )) %>%
#   mutate(pct = (pct.1 + pct.2) / 2) %>% 
#   ggplot(aes(pct, avg_log2FC)) +
#   geom_point(aes(color = sig), alpha = .8) +
#   geom_label_repel(aes(label = lbl), min.segment.length = .1) +
#   theme_bw() +
#   labs(x = "% Cells Expressing Gene", y = "Log2 Fold Change") +
#   facet_wrap(~subgroup, scales = "free_y") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   guides(color = guide_legend(""))
# 
# basal_genes <- de %>%
#   filter(subgroup == "Basal") %>%
#   pull(id)
# cycling_genes <- de %>%
#   filter(subgroup == "Cycling") %>%
#   pull(id)
# diff_genes <- de %>%
#   filter(subgroup == "Differentiating") %>%
#   pull(id)
# tsk_genes <- de %>%
#   filter(subgroup == "TSK") %>%
#   pull(id)
# genelist <- list("Basal" = basal_genes,
#                  "Cycling" = cycling_genes,
#                  "Differentiating" = diff_genes,
#                  "TSK" = tsk_genes)
