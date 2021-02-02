library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

data_dir <- file.path("data", "ji")
multiple_reps <- c("P1_cSCC", "P1_normal", "P3_cSCC", "P8_cSCC", "P8_normal")

samples <- list.dirs(data_dir, recursive=F)
# samples <- samples[1:5]

sample_list <- list()
for (i in 1:length(multiple_reps)) {
  s1name <- str_c(multiple_reps[i], "1", sep="_")
  s2name <- str_c(multiple_reps[i], "2", sep="_")
  s1_fp <- file.path(data_dir, s1name, "outs", "filtered_feature_bc_matrix")
  s2_fp <- file.path(data_dir, s2name, "outs", "filtered_feature_bc_matrix")
  s1 <- CreateSeuratObject(
    counts = Read10X(data.dir = s1_fp)
  )
  s2 <- CreateSeuratObject(
    counts = Read10X(data.dir = s2_fp)
  )
  print(s1)
  print(s2)
  s <- merge(s1, y = s2, project = multiple_reps[i])
  if (str_detect(s@project.name, "normal")) {
    s <- AddMetaData(s, "normal", col.name = "condition")
  } else {
    s <- AddMetaData(s, "cSCC", col.name = "condition")
  }
  sample_list[[i]] <- s
}

ji <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)])
rm(sample_list, s1, s2, s)

pattern <- str_c(multiple_reps, collapse="|")
for (i in 1:length(samples)) {
  if (!str_detect(samples[i], pattern)) {
    fp <- file.path(samples[i], "outs", "filtered_feature_bc_matrix")
    sname <- basename(samples[i])
    sname <- str_sub(sname, end = length(sname) - str_locate(rev(sname), "_")[1]+1)
    s <- CreateSeuratObject(
      counts = Read10X(data.dir = fp),
      project = sname
    )
    if (str_detect(s@project.name, "normal")) {
      s <- AddMetaData(s, "normal", col.name = "condition")
    } else {
      s <- AddMetaData(s, "cSCC", col.name = "condition")
    }
    ji <- merge(ji, y = s)
  }
}

rm(s, s1, s2)

# ji <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)])
# rm(sample_list)

ji[["percent.mt"]] <- PercentageFeatureSet(ji, pattern = "^MT-")
nfeat_plot <- VlnPlot(ji, "nFeature_RNA") + NoLegend()
ncount_plot <- VlnPlot(ji, "nCount_RNA") + NoLegend()
mt_plot <- VlnPlot(ji, "percent.mt",) + NoLegend()

print(nfeat_plot)
print(ncount_plot)
print(mt_plot)

ji <- subset(ji, nFeature_RNA > 200 & percent.mt < 10)

ji <- NormalizeData(ji)
ji <- FindVariableFeatures(ji, selection.method = "vst", nfeatures = 2000)
ji <- ScaleData(ji)

ji <- ji %>%
  RunPCA(features = VariableFeatures(ji)) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = .1) %>%
  RunUMAP(reduction = "pca", dims = 1:15)

DimPlot(ji, reduction = "umap")
