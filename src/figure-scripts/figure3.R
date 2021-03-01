# File: figure3.R
# Description: Generate plots for figure 2 of report

library(Seurat)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(fgsea)
library(patchwork)
library(ggrepel)
library(stringr)

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
p53_path <- gmtPathways("data/genesets/KEGG_P53_SIGNALING_PATHWAY.gmt")
emt_path <- gmtPathways("data/genesets/HALLMARK_EMT.gmt")
paths <- c(myc_path, p53_path, emt_path)

basal_de <- read_csv(basal_de_fp)
cycling_de <- read_csv(cycling_de_fp)
diff_de <- read_csv(diff_de_fp)
tsk_de <- read_csv(tsk_de_fp)

kcs <- AddModuleScore(kcs, features = paths, name = c("MYC_TARGETS", "P53_PATHWAY", "EMT_SCORE"))

kcs@meta.data %>%
  ggplot(aes(MYC_TARGETS1)) +
  geom_density(aes(fill = condition), alpha=.7) +
  theme_bw() +
  labs(x = "MYC TARGETS ACTIVITY", y = "Density") +
  scale_fill_discrete(name = "", labels=c("Non-PNI", "PNI")) +
  facet_wrap(~celltype)

kcs@meta.data %>%
  ggplot(aes(P53_PATHWAY2)) +
  geom_density(aes(fill = condition), alpha=.7) +
  theme_bw() +
  labs(x = "P53 PATHWAY ACTIVITY", y = "Density") +
  scale_fill_discrete(name = "", labels=c("Non-PNI", "PNI")) +
  facet_wrap(~celltype)

myc_plot <- kcs@meta.data %>%
  filter(celltype == "KC_Diff") %>%
  mutate(celltype = "Differentiating KCs") %>%
  ggplot(aes(MYC_TARGETS1)) +
  geom_density(aes(fill = condition), alpha=.7) +
  theme_bw() +
  labs(x = "MYC TARGETS ACTIVITY", y = "Density") +
  scale_fill_discrete(name = "", labels=c("Non-PNI", "PNI")) +
  facet_wrap(~celltype)

p53_plot <- kcs@meta.data %>%
  filter(celltype == "KC_Diff") %>%
  ggplot(aes(P53_PATHWAY2)) +
  geom_density(aes(fill = condition), alpha=.7) +
  theme_bw() +
  labs(x = "P53 PATHWAY ACTIVITY", y = "Density") +
  scale_fill_discrete(name = "", labels=c("Non-PNI", "PNI"))

activity_plots <- myc_plot / p53_plot + plot_layout(guides = "collect")
print(activity_plots)

## DE Analysis of KC Subtypes
de <- basal_de %>%
  bind_rows(cycling_de) %>%
  bind_rows(diff_de) %>%
  bind_rows(tsk_de) %>%
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) %>%
  mutate(id = str_c(gene, direction, sep = "_"))

de %>%
  mutate(sig = ifelse(p_val_adj < .05, "padj < .05", "ns")) %>%
  mutate(sig = factor(sig, levels = c("padj < .05", "ns"))) %>%
  mutate(lbl = case_when(
    subgroup == "Basal" & (avg_log2FC > 2 | avg_log2FC < -2) & p_val_adj < .05 ~ gene,
    subgroup == "Cycling" & (avg_log2FC > 2 | avg_log2FC < -2) & p_val_adj < .05 ~ gene,
    subgroup == "Differentiating" & (avg_log2FC > 4 | avg_log2FC < -2) & p_val_adj < .05 ~ gene,
    subgroup == "TSK" & (avg_log2FC > 2.5 | avg_log2FC < -2.3) & p_val_adj < .05 ~ gene,
    TRUE ~ ""
  )) %>%
  mutate(pct = (pct.1 + pct.2) / 2) %>% 
  ggplot(aes(pct, avg_log2FC)) +
  geom_point(aes(color = sig), alpha = .8) +
  geom_label_repel(aes(label = lbl), min.segment.length = .1) +
  theme_bw() +
  labs(x = "% Cells Expressing Gene", y = "Log2 Fold Change") +
  facet_wrap(~subgroup, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  guides(color = guide_legend(""))

basal_genes <- de %>%
  filter(subgroup == "Basal") %>%
  pull(id)
cycling_genes <- de %>%
  filter(subgroup == "Cycling") %>%
  pull(id)
diff_genes <- de %>%
  filter(subgroup == "Differentiating") %>%
  pull(id)
tsk_genes <- de %>%
  filter(subgroup == "TSK") %>%
  pull(id)
genelist <- list("Basal" = basal_genes,
                 "Cycling" = cycling_genes,
                 "Differentiating" = diff_genes,
                 "TSK" = tsk_genes)
