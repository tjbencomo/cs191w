# File: cellchat_kc_analysis.R
# Description: Analyze cellchat results when run with KC subpopulations
# in PNI data

library(CellChat)

infp <- file.path("data", "seurat", "harmony", "cellchat_pni_kc_specific.rds")
cellchat <- readRDS(infp)