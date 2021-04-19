## Consolidate metadata from broad labels and fine grained
## clustering into a single metadata file. 


library(Seurat)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble, include.only = "rownames_to_column")

data_dir <- file.path("data", "seurat", "harmony-20")
epiMetFp <- file.path(data_dir, "epithelial_metadata.csv")
endoMetFp <- file.path(data_dir, "endothelial_metadata.csv")
fibroMetFp <- file.path(data_dir, "fibroblast_metadata.csv")
melanMetFp <- file.path(data_dir, "melanocyte_metadata.csv")
myeloidMetFp <- file.path(data_dir, "myeloid_metadata.csv")
tcellMetFp <- file.path(data_dir, "tcell_metadata.csv")
cellsFp <- file.path(data_dir, "harmony_combined.rds")
outfp <- file.path(data_dir, "full_metadata.csv")

cells <- readRDS(cellsFp)
epiMet <- read_csv(epiMetFp) %>%
  select(-starts_with("RNA"), -seurat_clusters)
endoMet <- read_csv(endoMetFp) %>%
  select(-starts_with("RNA"), -seurat_clusters)
fibroMet <- read_csv(fibroMetFp) %>%
  select(-starts_with("RNA"), -seurat_clusters)
melanMet <- read_csv(melanMetFp) %>%
  select(-starts_with("RNA"), -seurat_clusters)
myeloidMet <- read_csv(myeloidMetFp) %>%
  select(-starts_with("RNA"), -seurat_clusters)
tcellMet <- read_csv(tcellMetFp) %>%
  select(-starts_with("RNA"), -seurat_clusters)

broadMet <- cells@meta.data %>%
  rownames_to_column(var = "cell_id") %>%
  as_tibble() %>%
  mutate(Broad_Label = Idents(cells)) %>%
  select(-starts_with("RNA"), -seurat_clusters)

join_cols <- c("cell_id", "orig.ident", "nCount_RNA", "nFeature_RNA", 
               "percent.mt", "condition")

met <- broadMet %>%
  left_join(epiMet, by = join_cols) %>%
  left_join(endoMet, by = join_cols) %>%
  left_join(fibroMet, by = join_cols) %>%
  left_join(melanMet, by = join_cols) %>%
  left_join(myeloidMet, by = join_cols) %>%
  left_join(tcellMet, by = join_cols)
# get all columns with keep information
keep_cols <- colnames(met)
keep_cols <- keep_cols[str_detect(keep_cols, "keep")]
subpop_cols <- colnames(met)
subpop_cols <- subpop_cols[str_detect(subpop_cols, "Label") & !subpop_cols == "Broad_Label"]

## Combine keep columns from subpop metadata
## Also combine subpopulation labels into a single column
met <- met %>%
  unite("keep_final", all_of(keep_cols), remove = T, na.rm = T) %>%
  unite("subpop_label", all_of(subpop_cols), remove = T, na.rm = T) %>%
  rename(
    Subpopulation_Label = subpop_label,
    keep = keep_final
  )

## B & Plasma cells + Mast cells were not subpop clustered
## Keep all of them
met <- met %>%
  mutate(keep = case_when(
    Broad_Label %in% c("Mast cells", "B & Plasma cells") ~ "Yes",
    TRUE ~ keep
  ))

## Add subpopulation label for B/Plasma and Mast cells - same as broad label
## Change broad label (or add a new column) for B & Plasma cells that were missclassifed
## as myeloid cells
## NOTE: Cell_Type is the finalized broad category cell type label
## Broad_Label is the label from the original broad UMAP clustering, which includes
## misclassifications like some B/Plasma cells called as Myeloid

met <- met %>%
  mutate(Subpopulation_Label = case_when(
    Broad_Label == "Mast cells" ~ "Mast cells",
    Broad_Label == "B & Plasma cells" ~ "B & Plasma cells",
    TRUE ~ Subpopulation_Label
  )) %>%
  mutate(Broad_Label = as.character(Broad_Label)) %>%
  mutate(Cell_Type = case_when(
    Broad_Label == "Myeloid cells" & Subpopulation_Label == "B & Plasma cells" ~ "B & Plasma cells",
    TRUE ~ Broad_Label
  ))

write_csv(met, outfp)  
