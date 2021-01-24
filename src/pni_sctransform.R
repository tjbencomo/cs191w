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
# rm(p1, p2, p3, p4)
# pni.list <- list(p1, p2, p3, p4)
pni.list <- list(p1, p4)

for (i in 1:length(pni.list)) {
  pni.list[[i]][["percent.mt"]] <- PercentageFeatureSet(pni.list[[i]], pattern = "^MT-")
  pni.list[[i]] <- subset(pni.list[[i]], nFeature_RNA > 200 & percent.mt < 20)
  pni.list[[i]] <- SCTransform(pni.list[[i]])
}

options(future.globals.maxSize= 1200*1024^2)
pni.features <- SelectIntegrationFeatures(object.list = pni.list, nfeatures = 3000)
pni.list <- PrepSCTIntegration(object.list = pni.list, anchor.features = pni.features, 
                               verbose = FALSE)


pni.anchors <- FindIntegrationAnchors(object.list = pni.list, normalization.method = "SCT", 
                       anchor.features = pni.list)
pni.integrated <- IntegrateData(anchorset = pni.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

