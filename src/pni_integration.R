library(Seurat)
library(patchwork)

data_dir <- "data"

patients <- list.files(data_dir)

p1 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE01", "filtered_feature_bc_matrix")),
  project = "LEE01"
)

p2 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE02", "filtered_feature_bc_matrix")),
  project = "LEE02"
)
p3 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE03", "filtered_feature_bc_matrix")),
  project = "LEE03"
)
p4 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE04", "filtered_feature_bc_matrix")),
  project = "LEE04"
)

batch2 <- merge(p2, y = list(p3, p4), add.cell.ids = patients[2:4])

p1 <- AddMetaData(p1, metadata = "batch1", col.name = "batch")
batch2 <- AddMetaData(batch2, metadata = "batch2", col.name = "batch")

p1[["percent.mt"]] <- PercentageFeatureSet(p1, pattern = "^MT-")
batch2[["percent.mt"]] <- PercentageFeatureSet(batch2, pattern = "^MT-")

p1 <- subset(p1, subset = nFeature_RNA > 200 & percent.mt < 20)
batch2 <- subset(batch2, subset = nFeature_RNA > 200 & percent.mt < 20)

p1 <- NormalizeData(p1)
p1 <- FindVariableFeatures(p1, selection.method = "vst", nfeatures = 2000)
batch2 <- NormalizeData(batch2)
batch2 <- FindVariableFeatures(batch2, selection.method = "vst", nfeatures = 2000)

pni.anchors <- FindIntegrationAnchors(object.list = list(p1, batch2), dims = 1:30)

pni <- IntegrateData(anchorset = pni.anchors, dims = 1:30)

DefaultAssay(pni) <- "integrated"
pni <- ScaleData(pni)
pni <- RunPCA(pni, features = VariableFeatures(object = pni))
pni <- FindNeighbors(pni, dims = 1:20)
pni <- FindClusters(pni, resolution = 0.1)
pni <- RunUMAP(pni, dims = 1:20)

cluster_plot <- DimPlot(pni, reduction = "umap")
patient_plot <- DimPlot(pni, reduction = "umap", group.by = "orig.ident")

seurat_integration <- cluster_plot + patient_plot +
  plot_annotation(title = "Seurat Integration")

## Try clustering without integration

cluster_no_integration <- function() {
  p1 <- CreateSeuratObject(
    counts = Read10X(data.dir = file.path(data_dir, "LEE01", "filtered_feature_bc_matrix")),
    project = "LEE01"
  )
  p2 <- CreateSeuratObject(
    counts = Read10X(data.dir = file.path(data_dir, "LEE02", "filtered_feature_bc_matrix")),
    project = "LEE02"
  )
  p3 <- CreateSeuratObject(
    counts = Read10X(data.dir = file.path(data_dir, "LEE03", "filtered_feature_bc_matrix")),
    project = "LEE03"
  )
  p4 <- CreateSeuratObject(
    counts = Read10X(data.dir = file.path(data_dir, "LEE04", "filtered_feature_bc_matrix")),
    project = "LEE04"
  )
  
  pni <- merge(p1, y = list(p2, p3, p4), add.cell.ids = patients, project = "PNI-SCC")
  rm(p1, p2, p3, p4)
  batch_info <- c(rep("batch1", sum(pni$orig.ident == "LEE01")), rep("batch2", sum(pni$orig.ident != "LEE01")))
  pni <- AddMetaData(pni, metadata = batch_info, col.name = "batch")
  
  pni[["percent.mt"]] <- PercentageFeatureSet(pni, pattern = "^MT-")
  
  pni <- subset(pni, subset = nFeature_RNA > 200 & percent.mt < 15)
  
  pni <- NormalizeData(pni)
  pni <- FindVariableFeatures(pni, selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(pni)
  pni <- ScaleData(pni, features = all.genes)
  pni <- RunPCA(pni, features = VariableFeatures(object = pni))
  
  pni <- FindNeighbors(pni, dims = 1:20)
  pni <- FindClusters(pni, resolution = 0.1)
  pni <- RunUMAP(pni, dims = 1:20)
  cluster_plot <- DimPlot(pni, reduction = "umap")
  patient_plot <- DimPlot(pni, reduction = "umap", group.by = "orig.ident")
  cluster_plot + patient_plot + plot_annotation(title = "No Integration")
}

no_int_plots <- cluster_no_integration()

seurat_integration / no_int_plots
