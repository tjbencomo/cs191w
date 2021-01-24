library(Seurat)

data_dir <- "data"

# p2 <- CreateSeuratObject(
#   counts = Read10X(data.dir = file.path(data_dir, "LEE02", "filtered_feature_bc_matrix")),
#   project = "LEE02"
# )
# 
# p2[["percent.mt"]] <- PercentageFeatureSet(p2, pattern = "^MT-")
# p2 <- subset(p2, nFeature_RNA > 200 & percent.mt < 20)
# p2 <- SCTransform(p2)


p1 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE01", "filtered_feature_bc_matrix")),
  project = "LEE01"
)
p3 <- CreateSeuratObject(
  counts = Read10X(data.dir = file.path(data_dir, "LEE03", "filtered_feature_bc_matrix")),
  project = "LEE03"
)

pni.list <- list(p1, p3)
rm(p1, p3)

for (i in 1:length(pni.list)) {
  pni.list[[i]][["percent.mt"]] <- PercentageFeatureSet(pni.list[[i]], pattern = "^MT-")
  pni.list[[i]] <- subset(pni.list[[i]], nFeature_RNA > 200 & percent.mt < 20)
  pni.list[[i]] <- SCTransform(pni.list[[i]])
}

pni.features <- SelectIntegrationFeatures(object.list = pni.list, nfeatures = 3000)
pni.list <- PrepSCTIntegration(object.list = pni.list, anchor.features = pni.features)

pni.anchors <- FindIntegrationAnchors(object.list = pni.list, normalization.method = "SCT", 
                                      anchor.features = pni.list)
pni.integrated <- IntegrateData(anchorset = pni.anchors, normalization.method = "SCT")