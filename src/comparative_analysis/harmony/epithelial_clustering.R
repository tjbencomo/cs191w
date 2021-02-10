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

data_dir <- "data"
fp <- file.path(data_dir, "seurat", "harmony", "harmony_combined.rds")
cells <- readRDS(fp)

epi <- subset(cells, idents = "Epithelial cells")
rm(cells)
print("Aftering excluding non-epithelial cells:")
print(table(epi$orig.ident))

epi <- epi %>%
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
non_kc_clusters <- c(1, 8, 15, 17)
kcs <- subset(epi, idents = non_kc_clusters, invert = T)
kc_outfp <- file.path(data_dir, "seurat", "harmony", "keratinocytes.rds")
saveRDS(kcs, kc_outfp)

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
  RunPCA(npc = 20) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony")
pni <- FindClusters(pni, resolution = .5)

pni_umap <- DimPlot(pni, reduction = "umap", label = T, label.size = 8)
pni_samples <- DimPlot(pni, reduction = "umap", group.by = "orig.ident")
pni_plots <- pni_umap + pni_samples

nonpni <- nonpni %>%
  RunPCA(npc = 20) %>%
  RunHarmony("orig.ident", max.iter.harmony = 30) %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony")
nonpni <- FindClusters(nonpni, resolution = .25)

nonpni_umap <- DimPlot(nonpni, reduction = "umap", label = T)
nonpni_samples <- DimPlot(nonpni, reduction = "umap", group.by = "orig.ident")
nonpni_plots <- nonpni_umap + nonpni_samples


(pni_umap + pni_samples) / (nonpni_umap + nonpni_samples)

pni_markers <- FindAllMarkers(pni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5) %>%
  filter(p_val_adj < .05)
top_pni_markers <- pni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()

nonpni_markers <- FindAllMarkers(nonpni, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5) %>%
  filter(p_val_adj < .05)
top_nonpni_markers <- nonpni_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  ungroup()


check_enrichr <- function(genes) {
  require(enrichR)
  setEnrichrSite("Enrichr")
  dbs <- c("MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", 
           "KEGG_2019_Human", "Reactome_2016", "GO_Molecular_Function_2018", 
           "GO_Biological_Process_2018")
  websiteLive <- TRUE
  if (is.null(dbs)) websiteLive <- FALSE
  if (websiteLive) {
    enriched <- enrichr(genes, dbs)
    bind_rows(enriched, .id = "db")
  } else {
    print("enrichR website not live!")
    stop()
  }
}

cluster1_res <- check_enrichr(pni_markers %>% filter(cluster == 1) %>% pull(gene))
cluster6_res <- check_enrichr(pni_markers %>% filter(cluster == 6) %>% pull(gene))
cluster2_res <- check_enrichr(pni_markers %>% filter(cluster == 2) %>% pull(gene))
cscc_basal_res <- check_enrichr(nonpni_markers %>% filter(cluster == 1) %>% pull(gene))

# kcs <- kcs %>%
#   RunPCA(npcs = 20) %>%
#   RunHarmony("condition", max.iter.harmony = 30) %>%
#   FindNeighbors(reduction = "harmony") %>%
#   FindClusters(resolution = .25) %>%
#   RunUMAP(reduction = "harmony", dims = 1:20)
# DimPlot(kcs, reduction = "umap", label = T)
# DimPlot(kcs, reduction = "umap", group.by = "condition")
# DimPlot(kcs, reduction = "umap", group.by = "orig.ident")
# # rm(epi)
# 
# kc.markers <- FindAllMarkers(kcs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = .5)
# kc.markers <- kc.markers %>%
#   filter(p_val_adj < .05)
# top_markers <- kc.markers %>%
#   group_by(cluster) %>%
#   slice_max(avg_log2FC, n = 15) %>%
#   ungroup()


# tsks <- subset(kcs, idents = 6)
# Idents(tsks) <- "condition"
# 
# avg.tsk.cells <- as.data.frame(log1p(AverageExpression(tsks, verbose = FALSE)$RNA))
# avg.tsk.cells$gene <- rownames(avg.tsk.cells)
# avg.tsk.cells$diff <- avg.tsk.cells$pni - avg.tsk.cells$cSCC
# 
# 
# degs <- FindMarkers(
#   tsks, 
#   ident.1 = "pni", 
#   ident.2 = "cSCC", 
#   verbose = TRUE, 
#   group.by="condition",
#   latent.vars = c("orig.ident", "percent.mt"),
#   test.use = "MAST"
#   ) %>%
#   tibble::rownames_to_column(var = "gene") %>%
#   filter(p_val_adj < .01) %>%
#   arrange(desc(abs(avg_log2FC)))
# 
# genes_to_label <- head(degs$gene, n=20)
# deg.plot <- avg.tsk.cells %>%
#   ggplot(aes(cSCC, pni)) +
#   geom_point() +
#   theme_bw() +
#   labs("cSCCs WITHOUT PNI", "cSCCs with PNI")
# deg.plot <- LabelPoints(plot = deg.plot, points = genes_to_label, repel = TRUE)
# deg.plot



