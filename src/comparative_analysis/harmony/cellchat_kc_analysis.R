# File: cellchat_kc_analysis.R
# Description: Analyze cellchat results when run with KC subpopulations
# in PNI data

library(CellChat)
library(stringr)

infp <- file.path("data", "seurat", "harmony", "cellchat_pni_kc_specific.rds")
cellchat <- readRDS(infp)

df.net <- subsetCommunication(cellchat)
df.net.pathways <- subsetCommunication(cellchat, slot.name = "netP")



groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

## See if any pathways unique to KC_Novel
df.net.pathways$id <- str_c(df.net.pathways$target, df.net.pathways$pathway_name, sep = ":")
kc_novel_paths <- df.net.pathways %>% filter(source == "KC_Novel") %>% pull(id)
other_kc_paths <- df.net.pathways %>% filter(str_detect(source, "KC") & !str_detect(source, "Novel")) %>% pull(id)
other_paths <- df.net.pathways %>% filter(!str_detect(source, "KC")) %>% pull(id)
print(setdiff(kc_novel_paths, other_kc_paths))
print(setdiff(kc_novel_paths, other_paths))


plot_pathway <- function(cellchat, pathway) {
  pathways.show <- c(pathway) 
  par(mfrow=c(1,1))
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(4,8) # a numeric vector. 
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
}

plot_pathway(cellchat, "NRG")
plot_pathway(cellchat, "GDNF")
plot_pathway(cellchat, "NGF")

## 
library(NMF)
library(ggalluvial)

# selectK(cellchat, pattern = "outgoing")

nPatternsOut = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatternsOut)

netAnalysis_river(cellchat, pattern = "outgoing")

# selectK(cellchat, pattern = "incoming")

nPatternsIn = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatternsIn)

netAnalysis_river(cellchat, pattern = "incoming")



