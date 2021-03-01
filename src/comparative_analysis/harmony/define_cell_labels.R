# File: define_cell_labels.R
# Description: Create metadata file with broad and fine-grained
# labels for each cell type. This can be used for CellPhoneDB
# analysis

library(Seurat)
library(dplyr)
library(readr)
library(stringr)
library(tibble)


broad_metadata_fp <- file.path("data", "combined_metadata.csv")
kcs_fp <- file.path("data", "seurat", "harmony", "keratinocytes.rds")
myeloid_fp <- file.path("data", "seurat", "harmony", "myeloid_cells.rds")
tcell_fp <- file.path("data", "seurat", "harmony", "T_cells.rds")
broad_metadata <- read_csv(broad_metadata_fp)

keep_cols <- c("cell_id", "orig.ident", "nCount_RNA", "nFeature_RNA", 
               "percent.mt", "condition", "seurat_clusters", "cell_type")

kcs <- readRDS(kcs_fp)
kcs@meta.data$cell_type <- Idents(kcs)
kc_metadata <- kcs@meta.data
kc_metadata <- kc_metadata %>%
  rownames_to_column(var = "cell_id") %>%
  select(all_of(keep_cols))
rm(kcs)

myeloid <- readRDS(myeloid_fp)
myeloid@meta.data$cell_type <- Idents(myeloid)
myeloid_metadata <- myeloid@meta.data
myeloid_metadata <- myeloid_metadata %>%
  rownames_to_column(var = "cell_id") %>%
  select(all_of(keep_cols))
rm(myeloid)

tcells <- readRDS(tcell_fp)
tcells@meta.data$cell_type <- Idents(tcells)
tcell_metadata <- tcells@meta.data
tcell_metadata <- tcell_metadata %>%
  rownames_to_column(var = "cell_id") %>%
  select(all_of(keep_cols))
rm(tcells)

fine_labels <- kc_metadata %>%
  bind_rows(myeloid_metadata) %>%
  bind_rows(tcell_metadata)

join_cols <- c("cell_id", "orig.ident", "nCount_RNA", "nFeature_RNA", 
               "percent.mt", "condition")
dc_types <- c("CLEC9A DCs", "CD1C DCs")


metadata <- broad_metadata %>%
  left_join(fine_labels, by = join_cols) %>%
  rename(
    celltype_level1 = broad_label,
    celltype_level3 = cell_type
  ) %>%
  mutate(celltype_level2 = case_when(
    celltype_level1 == "Epithelial cells" & is.na(celltype_level3) ~ "Unknown",
    celltype_level1 == "Epithelial cells" & !is.na(celltype_level3) ~ "Keratinocytes",
    celltype_level1 == "T cells/NK cells" & celltype_level3 == "NK" ~ "NK cells",
    celltype_level1 %in% c("T cells/NK cells", "Myeloid cells") & celltype_level3 == "Doublets" ~ "Doublets",
    celltype_level1 == "Myeloid cells" & celltype_level3 %in% dc_types ~ "DCs",
    str_detect(celltype_level3, "CD4") ~ "CD4+ T cells",
    str_detect(celltype_level3, "CD8") ~ "CD8+ T cells",
    str_detect(celltype_level3, "Treg") ~ "Treg",
    TRUE ~ celltype_level1
  )) %>%
    mutate(celltype_level3 = as.character(celltype_level3)) %>%
  mutate(celltype_level3 = case_when(
    !is.na(celltype_level3) ~ celltype_level3,
    TRUE ~ celltype_level2
  )) %>%
  mutate(isMultiplet = ifelse(celltype_level3 == "Doublets", 1, 0)) %>%
  select(-RNA_snn_res.0.6, -seurat_clusters.x, -seurat_clusters.y) %>%
  select(cell_id, orig.ident, nCount_RNA, nFeature_RNA, percent.mt, 
         condition, celltype_level1, celltype_level2, celltype_level3, 
         isMultiplet)


write_csv(metadata, "data/combined_cell_metadata.csv")    



# metadata <- broad_metadata %>%
#   left_join(fine_labels, by = join_cols) %>%
#   rename(celltype_level1 = broad_label,
#          celltype_level3 = cell_type) %>%
#   mutate(celltype_level3 = as.character(celltype_level3)) %>%
#   mutate(celltype_level3 = case_when(
#     is.na(celltype_level3) ~ "Not defined",
#     TRUE ~ celltype_level3
#   )) %>%
#   mutate(isMultiplet = ifelse(celltype_level3 == "Doublets", 1, 0)) %>%
#   mutate(celltype_level2 = case_when(
#     celltype_level1 == "Epithelial cells" & celltype_level3 == "Not defined" ~ "Pilosebaceous & Eccrine",
#     celltype_level1 == "Epithelial cells" & celltype_level3 != "Not defined" ~ "Keratinocyte",
#     celltype_level1 == "T cells/NK cells" & celltype_level3 == "NK" ~ "NK cells",
#     celltype_level1 == "T cells/NK cells" & celltype_level3 != "NK" & !isMultiplet ~ "T cells",
#     celltype_level1 == "Myeloid cells" & celltype_level3 %in% c("Macrophages", "Langerhans cells") ~ "Macrophages",
#     celltype_level1 == "Myeloid cells" & !(celltype_level3 %in% c("Macrophages", "Langerhans cells", "MDSCs")) & !isMultiplet ~ "DCs",
#     TRUE ~ celltype_level3
#   )) %>%
#   select(-RNA_snn_res.0.6, -seurat_clusters.x, -seurat_clusters.y) %>%
#   select(cell_id, orig.ident, nCount_RNA, nFeature_RNA, percent.mt, 
#          condition, celltype_level1, celltype_level2, celltype_level3, 
#          isMultiplet)

  
  
  
  
