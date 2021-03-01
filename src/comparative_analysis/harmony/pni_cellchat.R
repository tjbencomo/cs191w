# File: pni_cellchat.R
# Description: Perform crosstalk analysis on PNI single cell data

library(Seurat)
library(CellChat)

cells_fp <- file.path("data", "seurat", "harmony", "harmony_combined.rds")
metadata_fp <- file.path("data", "combined_cell_metadata.csv")

cells <- readRDS(cells_fp)
met <- readr::read_csv(metadata_fp)
lvl2_remove <- c("Unknown")
lvl3_remove <- c("KC_Unknown")
met <- met %>%
  dplyr::mutate(keep = ifelse(
    !celltype_level2 %in% lvl2_remove &
      !celltype_level3 %in% lvl3_remove &
      !(celltype_level1 == "Epithelial cells" & celltype_level3 == "Unknown") &
      isMultiplet != 1 &
      condition == "pni", 1, 0
  )) %>%
  dplyr::mutate(
    label = case_when(
      celltype_level2 == "Keratinocytes" ~ "KC",
      TRUE ~ celltype_level1
    )
  )
pni_met <- met %>%
  dplyr::filter(keep == 1)
cells$keep <- met$keep
pni <- subset(cells, keep == 1)
pni$cell_label <- pni_met$label
rm(cells)

cellchat <- createCellChat(object = pni, group.by = "cell_label")
cellChatDB <- CellChatDB.human
cellchat@DB <- cellChatDB
rm(pni)

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

## Visualize interactions
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1, 2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# mat <- cellchat@net$weight
# par(mfrow = c(3,2), xpd=TRUE)
# # for (i in 1:nrow(mat)) {
# for (i in 1:5) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }

# pathways.show <- c("MHC-I")
# vertex.receiver <- 4:8
# 
# netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
# 
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# 
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")


## Network analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

## Outgoing communication patterns
# library(NMF)
# library(ggalluvial)
# 
# selectK(cellchat, pattern = "outgoing")
# 
# nPatterns <- 2
# cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# 
# netAnalysis_river(cellchat, pattern = "outgoing")

print("Saving results object")
outfp <- file.path("data", "seurat", "harmony", "cellchat_pni.rds")
saveRDS(cellchat, outfp)


