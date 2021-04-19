# File: joint_cellchat.R
# Description: Compare and contrast cellular crosstalk in non-PNI and
# PNI samples


library(CellChat)
library(reticulate)
library(ComplexHeatmap)

use_python("/opt/anaconda3/envs/cellchat/bin/python")

data_dir <- file.path("data", "seurat", "harmony-20")
pni_fp <- file.path(data_dir, "cellchat_pni.rds")
nonpni_fp <- file.path(data_dir, "cellchat_ji.rds")

pni <- readRDS(pni_fp)
nonpni <- readRDS(nonpni_fp)

group.new = levels(pni@idents)
nonpni <- liftCellChat(nonpni, group.new)

object.list <- list("cSCC" = nonpni, "pSCC" = pni)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
rm(pni, nonpni)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) + 
  theme(text = element_text(size = 18))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") +
  theme(text = element_text(size = 14))
gg1 + gg2

netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

# group.cellType <- c("BC", "ENDO", "FIB", "KC", "MAST", "MELAN", "MYELOID", "T/NK")
# group.cellType <- factor(group.cellType, levels = group.cellType)
# object.list <- lapply(object.list, function(x) mergeInteractions(x, group.cellType))
# cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), 
                           attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)

rankSimilarity(cellchat, type = "functional")


# gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, cutoff.pvalue = 1e-3)
# gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
# gg1 + gg2
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = "KC_Novel-1", return.data = T) -> x
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, cutoff.pvalue = 1e-3)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, thresh = 1e-6, cutoff.pvalue = 1e-6)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, 
        sources.use = "KC")


df <- rankNet(cellchat, mode = "comparison", stacked = T, 
              do.stat = TRUE, sources.use = "KC", 
              return.data = T)[["signaling.contribution"]]
df %>% 
  filter(pvalues < 1e-4) %>% 
  ggplot(aes(name, contribution.scaled, fill=group)) + 
  geom_bar(stat = 'identity', position = 'fill') + 
  theme_classic() + 
  labs(x = "", y = "Relative Information Flow") +
  guides(fill = guide_legend(title = "")) +
  coord_flip() + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", size = 0.5) +
  theme(text = element_text(size = 18))




i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1,2),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(2),  comparison = c(1, 2), 
                        max.dataset = 2, title.name = "Increased signaling in PNI", angle.x = 45, 
                        remove.isolate = T, thresh = 1e-10, return.data = T)
gg1



pathways.show <- c("CLEC")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# weight.max <- getMaxWeight(object.list[[2]], slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))

## Visualize individual pathways
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle")
netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "circle")



