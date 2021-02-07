library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(harmony)

data_dir <- file.path("data", "ji")
data_dir <- "/scratch/users/tbencomo/cs191w/cellranger/ji"
multiple_reps <- c("P1_cSCC", "P1_normal", "P3_cSCC", "P8_cSCC", "P8_normal")

samples <- list.dirs(data_dir, recursive=F)
samples <- samples[str_detect(basename(samples), "P[0-9]+_cSCC")]
samples <- samples[!str_detect(basename(samples), "_2")]
print(samples)
# samples <- samples[1:5]

## Load Andrew's Samples
# sample_list <- list()
# for (i in 1:length(multiple_reps)) {
#   s1name <- str_c(multiple_reps[i], "1", sep="_")
#   s2name <- str_c(multiple_reps[i], "2", sep="_")
#   s1_fp <- file.path(data_dir, s1name, "outs", "filtered_feature_bc_matrix")
#   s2_fp <- file.path(data_dir, s2name, "outs", "filtered_feature_bc_matrix")
#   s1 <- CreateSeuratObject(
#     counts = Read10X(data.dir = s1_fp)
#   )
#   s2 <- CreateSeuratObject(
#     counts = Read10X(data.dir = s2_fp)
#   )
#   print(s1)
#   print(s2)
#   s <- merge(s1, y = s2, project = multiple_reps[i])
#   s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
#   s <- subset(s, nFeature_RNA > 200 & percent.mt < 10)
#   s <- SCTransform(s, vars.to.regress = "percent.mt", method = "glmGamPoi")
#   if (str_detect(s@project.name, "normal")) {
#     s <- AddMetaData(s, "normal", col.name = "condition")
#   } else {
#     s <- AddMetaData(s, "cSCC", col.name = "condition")
#   }
#   sample_list[[i]] <- s
# }

# ji <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)], merge.data = TRUE)
# rm(sample_list, s1, s2, s)

# print("Finished merging replicated samples")
# print("Moving on to remaining Ji samples...")


sample_list <- list()
pattern <- str_c(multiple_reps, collapse="|")
for (i in 1:length(samples)) {
  # if (!str_detect(samples[i], pattern)) {
    fp <- file.path(samples[i], "outs", "filtered_feature_bc_matrix")
    sname <- basename(samples[i])
    sname <- str_sub(sname, end = length(sname) - (str_locate(rev(sname), "_")[1])+1)
    s <- CreateSeuratObject(
      counts = Read10X(data.dir = fp),
      project = sname
    )
    s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
    s <- subset(s, nFeature_RNA > 200 & percent.mt < 10)
    # s <- NormalizeData(s)
    # s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
    # s <- SCTransform(s, vars.to.regress = "percent.mt", method = "glmGamPoi")
    if (str_detect(s@project.name, "normal")) {
      s <- AddMetaData(s, "normal", col.name = "condition")
    } else {
      s <- AddMetaData(s, "cSCC", col.name = "condition")
    }
    sample_list[[i]] <- s
    print(paste("Finished loading", sname))
    # ji <- merge(ji, y = s, merge.data = TRUE)
  # }
}

ji <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)], merge.data = TRUE)
# ji <- ji %>%
#     Normalize() %>%
#     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#     ScaleData(vars.to.regress = "percent.mt")
print("Finished loading and merging Ji samples")
print("Loading and merging PNI samples")

## Load PNI samples
pni_list <- list()
pni_dir <- "/scratch/users/tbencomo/cs191w/preprocessing/cellranger"
pni_samples <- list.dirs(pni_dir, recursive=F)

for (i in 1:length(pni_samples)) {
    sname <- basename(pni_samples[i])
    fp <- file.path(pni_samples[i], "outs", "filtered_feature_bc_matrix")
    p <- CreateSeuratObject(counts = Read10X(data.dir = fp), project = sname)
    p <- AddMetaData(p, "pni", col.name = "condition")
    p[["percent.mt"]] <- PercentageFeatureSet(p, pattern = "^MT-")
    p <- subset(p, nFeature_RNA > 200 & percent.mt < 20)
    # p <- SCTransform(p, vars.to.regress = "percent.mt", method = "glmGamPoi")
    pni_list[[i]] <- p
}

## Final Merge and Save
cells <- merge(ji, y = pni_list, merge.data = TRUE)
cells <- cells %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(vars.to.regress = "percent.mt") %>%
    RunPCA(npcs = 30)

cells <- RunHarmony(cells, "condition", max.iter.harmony = 30) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = .2)

markers <- FindAllMarkers(cells, only.pos = T, min.pct = .25, logfc.threshold = .5)
top_markers <- markers %>%
    filter(p_val_adj < .05) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10)

print(cells)
print(table(cells$orig.ident))

outdir = file.path("data", "seurat")
outfp <- file.path(outdir, "harmony_combined.rds")

saveRDS(cells, file = outfp)

write_csv(markers, file.path(outdir, "all_harmony_markers.csv"))
write_csv(top_markers, file.path(outdir, "top_harmony_markers.csv"))

