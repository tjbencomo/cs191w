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
# samples <- samples[1:5]
# print(samples)

## Load Andrew's Samples
multiple_reps <- multiple_reps[!str_detect(multiple_reps, "normal")]
rep_sample_list <- list()
for (i in 1:length(multiple_reps)) {
  s1name <- str_c(multiple_reps[i], "1", sep="_")
  s2name <- str_c(multiple_reps[i], "2", sep="_")
  s1_fp <- file.path(data_dir, s1name, "outs", "filtered_feature_bc_matrix")
  s2_fp <- file.path(data_dir, s2name, "outs", "filtered_feature_bc_matrix")
  s1 <- CreateSeuratObject(
    counts = Read10X(data.dir = s1_fp),
    project = multiple_reps[i]
  )
  s2 <- CreateSeuratObject(
    counts = Read10X(data.dir = s2_fp),
    project = multiple_reps[i]
  )
  # print(s1)
  # print(s2)
  s <- merge(s1, y = s2, project = multiple_reps[i])
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  s <- subset(s, nFeature_RNA > 200 & percent.mt < 10)
  if (str_detect(s@project.name, "normal")) {
    s <- AddMetaData(s, "normal", col.name = "condition")
  } else {
    s <- AddMetaData(s, "cSCC", col.name = "condition")
  }
  rep_sample_list[[i]] <- s
  print(s@project.name)
  print(table(s$orig.ident))
  print(s)
}

print("Samples with multiple replicates...")
print(rep_sample_list)


sample_list <- list()
idx <- 1
pattern <- str_c(multiple_reps, collapse="|")
for (i in 1:length(samples)) {
  if (!str_detect(samples[i], pattern)) {
    fp <- file.path(samples[i], "outs", "filtered_feature_bc_matrix")
    sname <- basename(samples[i])
    print(sname)
    sname <- str_c(str_split(sname, "_")[[1]][1:2], collapse = "_")
    # sname <- str_sub(sname, end = length(sname) - (str_locate(rev(sname), "_")[1])+1)
    s <- CreateSeuratObject(
      counts = Read10X(data.dir = fp),
      project = sname
    )
    s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
    s <- subset(s, nFeature_RNA > 200 & percent.mt < 10)
    if (str_detect(s@project.name, "normal")) {
      s <- AddMetaData(s, "normal", col.name = "condition")
    } else {
      s <- AddMetaData(s, "cSCC", col.name = "condition")
    }
    sample_list[[idx]] <- s
    print(paste("Finished loading", sname))
    print(sample_list[[idx]])
    idx <- idx + 1
  }
}

print("Merging sample lists")
# print(sample_list)
sample_list <- c(sample_list, rep_sample_list)
rm(rep_sample_list)
ji <- merge(sample_list[[1]], y = sample_list[2:length(sample_list)])


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
    pni_list[[i]] <- p
}

## Final Merge and Save
cells <- merge(ji, y = pni_list)
print(table(cells$orig.ident))

saveRDS(cells, "data/seurat/combined_raw.rds")
print("Saved Seurat object")

