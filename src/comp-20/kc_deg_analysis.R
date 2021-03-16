# Description: Perform differential expression analysis on KC subpopulations
# present in both tumor types

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(tibble, include.only = "rownames_to_column")

data_dir <- "data/seurat/harmony-20"
pni_fp <- file.path(data_dir, "pni_kcs.rds")
ji_fp <- file.path(data_dir, "ji_kcs.rds")

pni <- readRDS(pni_fp)
ji <- readRDS(ji_fp)

ji_labels <- c(
  `0` = "KC_Diff",
  `3` = "KC_Diff",
  `8` = "KC_Diff",
  `9` = "KC_Diff",
  `1` = "KC_Basal",
  `2` = "KC_Cycling",
  `4` = "KC_Cycling",
  `5` = "KC_TSK",
  `6` = "KC_U1",
  `7` = "KC_U2"
)
ji <- RenameIdents(ji, ji_labels)

pni_labels <- c(
  `0` = "KC_Basal",
  `2` = "KC_TSK",
  `5` = "KC_Cycling",
  `7` = "KC_Cycling",
  `9` = "KC_Cycling",
  `1` = "KC_Novel-1",
  `4` = "KC_Novel-2",
  `3` = "KC_Diff",
  `6` = "KC_Diff",
  `8` = "KC_Diff",
  `10` = "KC_U2"
)
pni <- RenameIdents(pni, pni_labels)

kcs <- merge(ji, pni, merge.data = T)

basal_de <- FindMarkers(
  kcs , ident.1 = "pni", ident.2 = "cSCC", 
  subset.ident = "KC_Basal", group.by = "condition",
  logfc.threshold = .5, min.pct = .25
  ) %>%
  tibble::rownames_to_column(var = "gene")

diff_de <- FindMarkers(
  kcs , ident.1 = "pni", ident.2 = "cSCC", 
  subset.ident = "KC_Diff", group.by = "condition",
  logfc.threshold = .5, min.pct = .25
 ) %>%
  tibble::rownames_to_column(var = "gene")

tsk_de <- FindMarkers(
  kcs , ident.1 = "pni", ident.2 = "cSCC", 
  subset.ident = "KC_TSK", group.by = "condition",
  logfc.threshold = .5, min.pct = .25
) %>%
  tibble::rownames_to_column(var = "gene")
