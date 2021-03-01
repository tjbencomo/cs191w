# File: epithelial_clustering.R
# Description: Isolate epithelial cells and identify
# subpopulations. We are looking for keratinocytes and
# eccrine (sweat gland) or pilosebaceous (hair/sebum) cells.
# We remove non-KCs for further KC analysis at the end.

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(ggrepel)

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

plotES <- function(es_res, db_fil = NULL, n_terms = 5) {
  s <- str_split(es_res$pathway, "_")
  db <- sapply(s, function(x) x[1])
  p <- sapply(s, function(x) str_c(x[1:length(x)], collapse =" "))
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


##############################################
## Analysis Code
##############################################

save_kcs <- FALSE
data_dir <- "data"
fp <- file.path(data_dir, "seurat", "harmony", "harmony_combined.rds")
cells <- readRDS(fp)

epi <- subset(cells, idents = "Epithelial cells")
rm(cells)
print("Aftering excluding non-epithelial cells:")
print(table(epi$orig.ident))

epi <- epi %>%
  NormalizeData() %>% # normally removed
  FindVariableFeatures() %>% # normally removed
  ScaleData(vars.to.regress = "percent.mt") %>% # normally removed
  RunPCA(npc = 30) %>%
  RunHarmony("condition", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony")
epi <- FindClusters(epi, resolution = .6)

p1 <- DimPlot(epi, reduction = "umap", label=T)
print(p1)
p2 <- DimPlot(epi, reduction = "umap", group.by = "condition")
print(p2)
p3 <- DimPlot(epi, reduction = "umap", group.by = "orig.ident")
print(p3)

pilosebaceous_featplot <- FeaturePlot(epi, c("SOX9", "LHX2", "KRT15", "KRT23", "SAA1"))
print(pilosebaceous_featplot)

epi.markers <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
epi.markers <- epi.markers %>%
  filter(p_val_adj < .05)
top_markers <- epi.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()
 
## We can now exclude pilosebaceous/eccrine clusters and proceed
## with KC analysis
# non_kc_clusters <- c(15, 17)
non_kc_clusters <- c(18)
kcs <- subset(epi, idents = non_kc_clusters, invert = T)

rm(epi)
pni <- subset(kcs, condition == "pni")
nonpni <- subset(kcs, condition == "cSCC")

basal_markers <- c("KRT15", "CCL2", "COL17A1", "CXCL14", 
                   "DST", "FTH1", "MT2A", "IGFBP5", "THBS2")
cycling_markers <- c("STMN1", "HIST1H4C", "TUBA1B", "PTTG1", 
                     "HMGB2", "H2AFZ", "TOP2A", "UBE2C", "NUSAP1")
diff_markers <- c("KRTDAP", "KRT10", "KRT1", "S100A7", "SBSN", 
                  "DMKN", "KRT60", "LYPD3", "KRT6A", "CALML5")
tsk_markers <- c("MMP10", "PTHLH", "FEZ1", "IL24", "KCNMA1",
                 "INHBA", "MAGEA4", "NT5E", "LAMC2", "SLITRK6")

pni <- pni %>%
  NormalizeData() %>% # normally removed
  FindVariableFeatures() %>% # normally removed
  ScaleData(vars.to.regress = "percent.mt") %>% # normally removed
  RunPCA(npc = 20) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony")
pni <- FindClusters(pni, resolution = .35)

pni_umap <- DimPlot(pni, reduction = "umap", label = T, label.size = 8)
print(pni_umap)
# pni_samples <- DimPlot(pni, reduction = "umap", group.by = "orig.ident")

pni_markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5) %>%
  filter(p_val_adj < .05)
top_pni_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()


## 8 = Remaining pilosebaceous/eccrine cluster
## 9 = Appear to be T cells that were improperly classified
## Remove pilosebaceous/t cells and recluster
pni_nonKC <- c(8, 9)
pni <- subset(pni, idents = pni_nonKC, invert = T)

pni <- pni %>%
  NormalizeData() %>% # normally removed
  FindVariableFeatures() %>% # normally removed
  ScaleData(vars.to.regress = "percent.mt") %>% # normally removed
  RunPCA(npc = 20) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony")
pni <- FindClusters(pni, resolution = .25)

pni_umap <- DimPlot(pni, reduction = "umap", label = T, label.size = 8)
print(pni_umap)

pni_markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5) %>%
  filter(p_val_adj < .05)
top_pni_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

pni_cluster_labels <- c(
  "KC_Diff",
  "KC_Novel",
  "KC_Cycling",
  "KC_Basal",
  "KC_Basal",
  "KC_TSK",
  "KC_Diff",
  "KC_Unknown"
)

# ## Cluster 3 is Unknown-1 and Cluster 5 is Unknown-2
# pni_cluster_labels <- c(
#   "Basal",
#   "Differentiating",
#   "Cycling",
#   "Unknown-1",
#   "TSK",
#   "Unknown-2",
#   "Differentiating"
# )

pni_cluster2id <- tibble(
  cluster = factor(0:(length(pni_cluster_labels)-1)),
  label = pni_cluster_labels
)

pni_markers <- pni_markers %>%
  left_join(pni_cluster2id)


names(pni_cluster_labels) <- levels(pni)
pni <- RenameIdents(pni, pni_cluster_labels)

nonpni <- nonpni %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 20) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony")
nonpni <- FindClusters(nonpni, resolution = .35)

nonpni_umap <- DimPlot(nonpni, reduction = "umap", label = T, label.size = 8)
print(nonpni_umap)
nonpni_samples <- DimPlot(nonpni, reduction = "umap", group.by = "orig.ident")
nonpni_plots <- nonpni_umap + nonpni_samples


nonpni_markers <- FindAllMarkers(nonpni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5) %>%
  filter(p_val_adj < .05)
top_nonpni_markers <- nonpni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

## 7 - Appear to be Myeloid KC doublets
## 10 - Outliers that don't cluster well and very few; also exclude
nonpni_nonKC <- c(7, 10)
nonpni <- subset(nonpni, idents = nonpni_nonKC, invert = T)

nonpni <- nonpni %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA(npc = 20) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony")
nonpni <- FindClusters(nonpni, resolution = .25)

nonpni_umap <- DimPlot(nonpni, reduction = "umap", label = T, label.size = 8)
print(nonpni_umap)

nonpni_markers <- FindAllMarkers(nonpni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5) %>%
  filter(p_val_adj < .05)
top_nonpni_markers <- nonpni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

# ## Cluster 6 is the unknown non-PNI KC cluster
# nonpni_cluster_labels <- c(
#   "Differentiating",
#   "Basal",
#   "Cycling",
#   "Differentiating",
#   "Cycling",
#   "TSK",
#   "Unknown-3"
# )

nonpni_cluster_labels <- c(
  "KC_Diff",
  "KC_Basal",
  "KC_Cycling",
  "KC_Diff",
  "KC_Cycling",
  "KC_TSK",
  "KC_Diff"
)

nonpni_cluster2id <- tibble(
  cluster = factor(0:(length(nonpni_cluster_labels)-1)),
  label = nonpni_cluster_labels
)

nonpni_markers <- nonpni_markers %>%
  left_join(nonpni_cluster2id)

names(nonpni_cluster_labels) <- levels(nonpni)
nonpni <- RenameIdents(nonpni, nonpni_cluster_labels)

## Redraw with cluster labels
pni_umap <- DimPlot(pni, reduction = "umap", label = T, label.size = 8)
nonpni_umap <- DimPlot(nonpni, reduction = "umap", label = T, label.size = 8)

(pni_umap + ggtitle("PNI Keratinocytes")) / (nonpni_umap + ggtitle("Non-PNI Keratinoctyes"))

## Re-merge cells now that we have labeled the clusters separately
kcs <- merge(nonpni, y = pni, merge.data = T)


novel.markers <- FindMarkers(pni, ident.1 = "KC_Novel", logfc.threshold = .25) %>% 
  tibble::rownames_to_column(var = "gene") %>%
  filter(p_val_adj < .05) %>%
  arrange(avg_log2FC )
# u2.markers <- FindMarkers(pni, ident.1 = "Unknown-2", logfc.threshold = .25) %>% 
#   tibble::rownames_to_column(var = "gene") %>%
#   filter(p_val_adj < .05) %>%
#   arrange(avg_log2FC)
# u3.markers <- FindMarkers(nonpni, ident.1 = "Unknown-3", logfc.threshold = .25) %>% 
#   tibble::rownames_to_column(var = "gene") %>%
#   filter(p_val_adj < .05) %>%
#   arrange(avg_log2FC)


novel_ranks <- novel.markers$avg_log2FC
names(novel_ranks) <- novel.markers$gene
novel_es <- run_fgsea(novel_ranks)
novel_esplot <- plotES(novel_es)
print(novel_esplot)

# u2_ranks <- u2.markers$avg_log2FC
# names(u2_ranks) <- u2.markers$gene
# u2_es <- run_fgsea(u2_ranks)
# 
# u3_ranks <- u3.markers$avg_log2FC
# names(u3_ranks) <- u3.markers$gene
# u3_es <- run_fgsea(u3_ranks)

# u1_esplot <- plotES(u1_es, n_terms = 6)
# print(u1_esplot)
# u2_esplot <- plotES(u2_es, n_terms = 6)
# print(u2_esplot)
# u3_esplot <- plotES(u3_es, n_terms = 6)
# print(u3_esplot)


kc_de <- FindMarkers(
  kcs,
  ident.1 = "pni",
  ident.2 = "cSCC",
  group.by = "condition",
  logfc.threshold = .5
  ) %>%
  tibble::rownames_to_column(var = "gene")

basal_de <- FindMarkers(
  kcs, ident.1 = "pni", 
  ident.2 = "cSCC",
  subset.ident = "KC_Basal",
  group.by = "condition", 
  logfc.threshold = .5) %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(subgroup = "Basal")

cycling_de <- FindMarkers(
  kcs, ident.1 = "pni", 
  ident.2 = "cSCC",
  subset.ident = "KC_Cycling",
  group.by = "condition", 
  logfc.threshold = .5) %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(subgroup = "Cycling")

diff_de <- FindMarkers(
  kcs, ident.1 = "pni", 
  ident.2 = "cSCC",
  subset.ident = "KC_Diff",
  group.by = "condition", 
  logfc.threshold = .5) %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(subgroup = "Differentiating")

tsk_de <- FindMarkers(
  kcs, ident.1 = "pni", 
  ident.2 = "cSCC",
  subset.ident = "KC_TSK",
  group.by = "condition", 
  logfc.threshold = .5) %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(subgroup = "TSK")


library(UpSetR)
de_upset <- upset(fromList(list(
  Basal = basal_de %>% filter(p_val_adj < .05) %>% pull(gene),
  Cycling = cycling_de %>% filter(p_val_adj < .05) %>% pull(gene),
  Differentiating = diff_de %>% filter(p_val_adj < .05) %>% pull(gene),
  TSK = tsk_de %>% filter(p_val_adj < .05) %>% pull(gene)
)), text.scale = 2, order.by = "freq")
print(de_upset)

de_plot <- bind_rows(basal_de, cycling_de, diff_de, tsk_de) %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(sig = ifelse(p_val_adj < .05, "padj < .05", "ns")) %>%
  ggplot(aes(avg_log2FC, logp)) +
  geom_point(aes(color = sig)) +
  facet_wrap(~subgroup, scales = "free")
print(de_plot)

basal_ranks <- basal_de %>%
  arrange(avg_log2FC) %>%
  pull(avg_log2FC)
names(basal_ranks) <- basal_de %>%
  arrange(avg_log2FC) %>%
  pull(gene)
basal_es <- run_fgsea(basal_ranks)
basal_esplot <- plotES(basal_es, n_terms = 4)
basal_de %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(sig = ifelse(p_val_adj < .05, "padj < .05", "ns")) %>%
  mutate(lbl = ifelse(abs(avg_log2FC) > 2, gene, "")) %>%
  ggplot(aes(avg_log2FC, logp, label = lbl)) +
  geom_point(aes(color = sig)) +
  geom_label_repel()

cycling_ranks <- cycling_de %>%
  arrange(avg_log2FC) %>%
  pull(avg_log2FC)
names(cycling_ranks) <- cycling_de %>%
  arrange(avg_log2FC) %>%
  pull(gene)
cycling_es <- run_fgsea(cycling_ranks)
cycling_esplot <- plotES(cycling_es, n_terms = 4)
cycling_de %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(sig = ifelse(p_val_adj < .05, "padj < .05", "ns")) %>%
  mutate(lbl = ifelse(abs(avg_log2FC) > 2, gene, "")) %>%
  ggplot(aes(avg_log2FC, logp, label = lbl)) +
  geom_point(aes(color = sig)) +
  geom_label_repel()


diff_ranks <- diff_de %>%
  arrange(avg_log2FC) %>%
  pull(avg_log2FC)
names(diff_ranks) <- diff_de %>%
  arrange(avg_log2FC) %>%
  pull(gene)
diff_es <- run_fgsea(diff_ranks)

diff_esplot <- plotES(diff_es, n_terms = 4)
diff_de %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(sig = ifelse(p_val_adj < .05, "padj < .05", "ns")) %>%
  mutate(lbl = ifelse(abs(avg_log2FC) > 3, gene, "")) %>%
  ggplot(aes(avg_log2FC, logp, label = lbl)) +
  geom_point(aes(color = sig)) +
  geom_label_repel()


tsk_ranks <- tsk_de %>%
  arrange(avg_log2FC) %>%
  pull(avg_log2FC)
names(tsk_ranks) <- tsk_de %>%
  arrange(avg_log2FC) %>%
  pull(gene)
tsk_es <- run_fgsea(tsk_ranks)
tsk_esplot <- plotES(tsk_es, n_terms = 6)
tsk_de %>%
  mutate(logp = -log10(p_val_adj)) %>%
  mutate(sig = ifelse(p_val_adj < .05, "padj < .05", "ns")) %>%
  mutate(lbl = ifelse(abs(avg_log2FC) > 2, gene, "")) %>%
  ggplot(aes(avg_log2FC, logp, label = lbl)) +
  geom_point(aes(color = sig)) +
  geom_label_repel()

## Save keratinocyte objects for future use
kc_outfp <- file.path(data_dir, "seurat", "harmony", "keratinocytes.rds")
pni_outfp <- file.path(data_dir, "seurat", "harmony", "pni_kcs.rds")
nonpni_outfp <- file.path(data_dir, "seurat", "harmony", "nonpni_kcs.rds")
if (save_kcs) {
  saveRDS(kcs, kc_outfp)
  saveRDS(pni, pni_outfp)
  saveRDS(nonpni, nonpni_outfp)
}

## Save markers genes for later viewing
write_csv(pni_markers, "data/seurat/harmony/pni_markers.csv")
write_csv(nonpni_markers, "data/seurat/harmony/nonpni_markers.csv")
write_csv(novel.markers, "data/seurat/harmony/novel_signature.csv")
write_csv(novel_es[, -"leadingEdge"], "data/seurat/harmony/novel_enrichment.csv")

## Save DE gene lists
write_csv(kc_de, "data/seurat/harmony/kc_de_res.csv")
write_csv(basal_de, "data/seurat/harmony/basal_de_res.csv")
write_csv(cycling_de, "data/seurat/harmony/cycling_de_res.csv")
write_csv(diff_de, "data/seurat/harmony/diff_de_res.csv")
write_csv(tsk_de, "data/seurat/harmony/tsk_de_res.csv")

