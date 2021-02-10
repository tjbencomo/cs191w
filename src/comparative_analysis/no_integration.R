library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(glmGamPoi)
library(purrr, include.only = "map_dfr")
library(tibble, include.only = "rownames_to_column")

data_dir <- file.path("data", "ji")
data_dir <- "/scratch/users/tbencomo/cs191w/cellranger/ji"
# multiple_reps <- c("P1_cSCC", "P1_normal", "P3_cSCC", "P8_cSCC", "P8_normal")

# samples <- list.dirs(data_dir, recursive=F)
# samples <- samples[str_detect(basename(samples), "P[0-9]+_cSCC")]
# samples <- samples[!str_detect(basename(samples), "_2")]
# print(samples)
# # samples <- samples[1:5]

# ## Load Andrew's Samples
# # sample_list <- list()
# # for (i in 1:length(multiple_reps)) {
# #   s1name <- str_c(multiple_reps[i], "1", sep="_")
# #   s2name <- str_c(multiple_reps[i], "2", sep="_")
# #   s1_fp <- file.path(data_dir, s1name, "outs", "filtered_feature_bc_matrix")
# #   s2_fp <- file.path(data_dir, s2name, "outs", "filtered_feature_bc_matrix")
# #   s1 <- CreateSeuratObject(
# #     counts = Read10X(data.dir = s1_fp)
# #   )
# #   s2 <- CreateSeuratObject(
# #     counts = Read10X(data.dir = s2_fp)
# #   )
# #   print(s1)
# #   print(s2)
# #   s <- merge(s1, y = s2, project = multiple_reps[i])
# #   s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
# #   s <- subset(s, nFeature_RNA > 200 & percent.mt < 10)
# #   s <- SCTransform(s, vars.to.regress = "percent.mt", method = "glmGamPoi")
# #   if (str_detect(s@project.name, "normal")) {
# #     s <- AddMetaData(s, "normal", col.name = "condition")
# #   } else {
# #     s <- AddMetaData(s, "cSCC", col.name = "condition")
# #   }
# #   sample_list[[i]] <- s
# # }

# # ji <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)], merge.data = TRUE)
# # rm(sample_list, s1, s2, s)

# print("Finished merging replicated samples")
# # print("Moving on to remaining Ji samples...")


# sample_list <- list()
# pattern <- str_c(multiple_reps, collapse="|")
# for (i in 1:length(samples)) {
#   # if (!str_detect(samples[i], pattern)) {
#     fp <- file.path(samples[i], "outs", "filtered_feature_bc_matrix")
#     sname <- basename(samples[i])
#     sname <- str_sub(sname, end = length(sname) - (str_locate(rev(sname), "_")[1])+1)
#     s <- CreateSeuratObject(
#       counts = Read10X(data.dir = fp),
#       project = sname
#     )
#     s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
#     s <- subset(s, nFeature_RNA > 200 & percent.mt < 10)
#     # s <- NormalizeData(s)
#     # s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
#     # s <- SCTransform(s, vars.to.regress = "percent.mt", method = "glmGamPoi")
#     if (str_detect(s@project.name, "normal")) {
#       s <- AddMetaData(s, "normal", col.name = "condition")
#     } else {
#       s <- AddMetaData(s, "cSCC", col.name = "condition")
#     }
#     sample_list[[i]] <- s
#     print(paste("Finished loading", sname))
#     # ji <- merge(ji, y = s, merge.data = TRUE)
#   # }
# }

# ji <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)], merge.data = TRUE)
# # ji <- ji %>%
# #     Normalize() %>%
# #     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
# #     ScaleData(vars.to.regress = "percent.mt")
# print("Finished loading and merging Ji samples")
# print("Loading and merging PNI samples")

# ## Load PNI samples
# pni_list <- list()
# pni_dir <- "/scratch/users/tbencomo/cs191w/preprocessing/cellranger"
# pni_samples <- list.dirs(pni_dir, recursive=F)

# for (i in 1:length(pni_samples)) {
#     sname <- basename(pni_samples[i])
#     fp <- file.path(pni_samples[i], "outs", "filtered_feature_bc_matrix")
#     p <- CreateSeuratObject(counts = Read10X(data.dir = fp), project = sname)
#     p <- AddMetaData(p, "pni", col.name = "condition")
#     p[["percent.mt"]] <- PercentageFeatureSet(p, pattern = "^MT-")
#     p <- subset(p, nFeature_RNA > 200 & percent.mt < 20)
#     # p <- SCTransform(p, vars.to.regress = "percent.mt", method = "glmGamPoi")
#     pni_list[[i]] <- p
# }

# ## Final Merge and Save
# cells <- merge(ji, y = pni_list, merge.data = TRUE)
infp <- file.path("data/seurat/combined_raw.rds")
cells <- readRDS(infp)
cells <- cells %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(vars.to.regress = "percent.mt") %>%
    RunPCA(npcs = 30) %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>%
    FindNeighbors() %>%
    FindClusters(resolution = .4)

p1 <- DimPlot(cells, reduction = "umap", label = T)
p2 <- DimPlot(cells, reduction = "umap", group.by = "condition")
p3 <- DimPlot(cells, reduction = "umap", group.by = "orig.ident")
combo_plot <- p1 + p2

ggsave("figures/noint_umap.eps", combo_plot, width = 12, height = 6)
ggsave("figures/noint_samples.eps", p3)




get_conserved <- function(cluster){
    FindConservedMarkers(cells,
    ident.1 = cluster,
    grouping.var = "condition",
    only.pos = TRUE, logfc.threshold = .5) %>%
    rownames_to_column(var = "gene") %>%
    mutate(cluster_id = cluster)
}

cluster_ids <- 0:(length(unique(cells$seurat_clusters))-1)
markers <- map_dfr(cluster_ids, get_conserved)
top_markers <- markers %>%
  mutate(combined_avg_log2FC = (pni_avg_log2FC + cSCC_avg_log2FC) / 2) %>%
  group_by(cluster_id) %>%
    slice_max(combined_avg_log2FC, n=10) %>%
    ungroup() %>%
    select(gene, combined_avg_log2FC, cluster_id, max_pval, everything())




print(cells)
print(table(cells$orig.ident))

outdir = file.path("data", "seurat")
outfp <- file.path(outdir, "no_int_combined.rds")

saveRDS(cells, file = outfp)

write_csv(markers, file.path(outdir, "all_noint_markers.csv"))
write_csv(top_markers, file.path(outdir, "top_noint_markers.csv"))

