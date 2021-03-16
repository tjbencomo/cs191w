# Description: Investigate the novel-1 cluster found in PNI samples

library(Seurat)
library(dplyr)
library(readr)
library(enrichR)
library(stringr)
library(data.table, include.only = "rbindlist")

data_dir <- file.path("data", "seurat", "harmony-20")
pni_marker_fp <- file.path(data_dir, "pni_kc_markers.csv")
pni_markers <- read_csv(pni_marker_fp)

markers <- pni_markers %>%
  filter(cluster == 1) %>%
  pull(gene)

dbs <- c("KEGG_2019_Human", "MSigDB_Hallmark_2020", 
         "WikiPathways_2019_Human", "GO_Biological_Process_2018", 
         "GO_Molecular_Function_2018", "Reactome_2016")
res <- enrichr(markers, dbs)
for (i in 1:length(res)) {
  res[[i]]$db <- names(res)[i]
}
resdf <- rbindlist(res)
resdf <- resdf %>%
  filter(Adjusted.P.value < .05)

hallmark <- resdf %>%
  filter(db == "MSigDB_Hallmark_2020") %>%
  slice_min(Adjusted.P.value, n = 9) %>%
  mutate(nGenes = as.numeric(str_split(Overlap, "/", simplify = T)[, 1]))

hallmark %>%
  mutate(cluster = "Novel-1") %>%
  ggplot(aes(cluster, reorder(Term, -Adjusted.P.value))) +
  geom_point(aes(color = -log10(Adjusted.P.value), size = nGenes)) +
  scale_size(range = c(2, 8)) +
  theme_bw() +
  labs(x = "", y = "")

kegg <- resdf %>%
  filter(db == "KEGG_2019_Human") %>%
  slice_min(Adjusted.P.value, n = 9) %>%
  mutate(nGenes = as.numeric(str_split(Overlap, "/", simplify = T)[, 1]))

kegg %>%
  tail(9) %>%
  mutate(cluster = "Novel-1") %>%
  ggplot(aes(cluster, reorder(Term, -Adjusted.P.value))) +
  geom_point(aes(color = -log10(Adjusted.P.value), size = nGenes)) +
  scale_size(range = c(2, 8)) +
  theme_bw() +
  labs(x = "", y = "")

wikipath <- resdf %>%
  filter(db == "WikiPathways_2019_Human") %>%
  slice_min(Adjusted.P.value, n = 9) %>%
  mutate(nGenes = as.numeric(str_split(Overlap, "/", simplify = T)[, 1]))

wikipath %>%
  tail(8) %>%
  mutate(cluster = "Novel-1") %>%
  ggplot(aes(cluster, reorder(Term, -Adjusted.P.value))) +
  geom_point(aes(color = -log10(Adjusted.P.value), size = nGenes)) +
  scale_size(range = c(2, 8)) +
  theme_bw() +
  labs(x = "", y = "")

pni_markers %>%
  filter(cluster == 1) %>%
  mutate(ribosomal = case_when(
    str_starts(gene, "RPL") ~ "Ribosomal",
    TRUE ~ "Other"
  )) %>%
  ggplot(aes(ribosomal)) +
  geom_bar(aes(fill = ribosomal)) +
  theme_bw() +
  labs(x = "", y = "Number of Genes") +
  theme(
    text = element_text(size = 14)
  ) +
  guides(fill = FALSE)
