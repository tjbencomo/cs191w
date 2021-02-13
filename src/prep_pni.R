library(Seurat)

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

pni <- merge(p1, y = list(p2, p3, p4), add.cell.ids = c("LEE01", "LEE02", "LEE03", "LEE04"))
rm(p1, p2, p3, p4)

pni[["percent.mt"]] <- PercentageFeatureSet(pni, pattern = "^MT-")

VlnPlot(pni, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

