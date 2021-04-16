# File: ji_cellchat.R
# Description: Perform crosstalk analysis on Ji's single cell data

library(Seurat)
library(readr)
library(dplyr)
library(CellChat)


data_dir <- file.path("data", "seurat", "harmony-20")
cells_fp <- file.path(data_dir, "harmony_combined.rds")
metadata_fp <- file.path(data_dir, "full_metadata.csv")
outfp <- file.path(data_dir, "cellchat_ji.rds")

cells <- readRDS(cells_fp)
metadata <- read_csv(metadata_fp)


## Collapse fibro/endo subcategories
cells2remove <- c("Pilosebaceous/Eccrine cell", "Melanocyte")
metadata <- metadata %>%
  filter(keep == "Yes", condition == "cSCC") %>%
  filter(!(Subpopulation_Label %in% cells2remove)) %>%
  mutate(Final_Label = case_when(
    Cell_Type == "Endothelial cells" ~ "Endo",
    Cell_Type == "Fibroblasts" ~ "Fibro",
    Cell_Type == "Mast cells" ~ "Mast",
    Cell_Type == "B & Plasma cells" ~ "B & Plasma",
    Cell_Type == "Myeloid cells" ~ "Myeloid",
    Cell_Type == "T/NK cells" ~ "T/NK",
    Cell_Type == "Epithelial cells" ~ Subpopulation_Label,
    TRUE ~ Cell_Type
  ))

cells$keep <- rownames(cells@meta.data) %in% metadata$cell_id
ji <- subset(cells, keep == 1)
rm(cells)
ji$cell_label <- metadata$Final_Label

cellchat <- createCellChat(object = ji, group.by = "cell_label")
cellChatDB <- CellChatDB.human
cellchat@DB <- cellChatDB
rm(ji)

if (Sys.getenv("RSTUDIO") != "1") {
  print("Using multiple cores")
  future::plan("multicore", workers = 3)
}


cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

print("Computing communication probabilities")
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
df.net.pathways <- subsetCommunication(cellchat, slot.name = "netP")

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

## Network Analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

print("Saving results object")
saveRDS(cellchat, outfp)
