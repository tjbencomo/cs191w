library(dplyr)
library(Seurat)
library(patchwork)

data_dir <- file.path("data", "LEE01", "filtered_feature_bc_matrix")
scc.data <- Read10X(data.dir = data_dir)

pniscc <- CreateSeuratObject(counts = scc.data, project = "scc1k", min.cells = 3, min.features = 200)
pniscc

pniscc[["percent.mt"]] <- PercentageFeatureSet(pniscc, pattern = "^MT-")

VlnPlot(pniscc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pniscc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pniscc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pniscc <- subset(pniscc, subset = nFeature_RNA > 200 & percent.mt < 15)
pniscc <- NormalizeData(pniscc)

pniscc <- FindVariableFeatures(pniscc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pniscc), 10)
plot1 <- VariableFeaturePlot(pniscc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pniscc)
pniscc <- ScaleData(pniscc, features = all.genes)
pniscc <- RunPCA(pniscc, features = VariableFeatures(object = pniscc))

VizDimLoadings(pniscc, dims = 1:2, reduction = "pca")
DimPlot(pniscc, reduction = "pca")

DimHeatmap(pniscc, dims = 1, cells = 500, balanced = TRUE)

pniscc <- JackStraw(pniscc, num.replicate = 100)
pniscc <- ScoreJackStraw(pniscc, dims = 1:20)
JackStrawPlot(pniscc, dims = 1:20)


ElbowPlot(pniscc, ndims = 50)


pniscc <- FindNeighbors(pniscc, dims = 1:30)
pniscc <- FindClusters(pniscc, resolution = 0.8)
pniscc <- RunUMAP(pniscc, dims = 1:15)
DimPlot(pniscc, reduction = "umap")

pniscc.markers <- FindAllMarkers(pniscc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pniscc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> markers


library(celldex)
library(SingleR)

pred.hesc <- SingleR(test = hESCs, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)
