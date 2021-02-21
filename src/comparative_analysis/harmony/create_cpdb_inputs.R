# File: create_cpdb_inputs.R
# Description: Create metadata and RNA abundance files for
# CellPhoneDB analysis. Analysis cSCCs and PNI samples separately
# for comparative analysis

library(Seurat)
library(readr)
library(dplyr)

cells_fp <- file.path("data", "seurat", "harmony", "harmony_combined.rds")
metadata_fp <- file.path("data", "combined_cell_metadata.csv")
cells <- readRDS(cells_fp)
metadata <- read_csv(metadata_fp)

isDoublet <- metadata$isMultiplet == 1
cells$isDoublet <- isDoublet

cells <- subset(cells, isDoublet == 0)

ji <- subset(cells, condition == "cSCC")
pni <- subset(cells, condition == "pni")
rm(cells)

ji_metadata <- metadata %>% 
  filter(condition == "cSCC", isMultiplet == 0)
pni_metadata <- metadata %>% 
  filter(condition == "pni", isMultiplet == 0)

print("Checking that cell IDs match metadata order...")
stopifnot(all(ji_metadata$cell_id == rownames(ji@meta.data)))
stopifnot(all(pni_metadata$cell_id == rownames(pni@meta.data)))
print("CellID check passed!")

print("Writing inputs for Ji cohort...")
ji_counts <- ji@assays$RNA@data %>%
  as_tibble() %>% 
  rownames_to_column(var = "Gene")
ji_celltypes <- ji_metadata %>%
  select(cell_id, celltype_level3)
write_csv(ji_counts, 'data/cpdb/ji_counts.txt')
write_tsv(ji_celltypes, 'data/cpdb/ji_meta.txt')

print("Writing inputs for PNI cohort")
pni_counts <- pni@assays$RNA@data %>%
  as_tibble() %>% 
  rownames_to_column(var = "Gene")
pni_celltypes <- pni_metadata %>%
  select(cell_id, celltype_level3)
write_csv(pni_counts, 'data/cpdb/pni_counts.csv')
write_csv(pni_celltypes, 'data/cpdb/pni_meta.csv')

