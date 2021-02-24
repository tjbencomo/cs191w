# File: compare_cell_barcodes.R
# Description: Investigate whether my preprocessing of Andrew's
# scRNA-seq data has retained more cells than his preprocessing run
# using an older version of cellranger

library(readr)
library(dplyr)
library(tibble)
library(stringr)

geo_cells_fp <- file.path("data", "GSE144236_patient_metadata_new.txt")
my_cells_fp <- file.path("data", "combined_metadata.csv")

geo_cells <- read_tsv(geo_cells_fp)
cells <- read_csv(my_cells_fp)

## Remove normal cells - we didn't look at those
geo_cells <- geo_cells %>%
  filter(tum.norm == "Tumor") %>%
  mutate(barcode = str_split(cell_id, "_|-", simplify = T)[, 3])


## Remove PNI cells
cells <- cells %>%
  filter(condition == "cSCC", nFeature_RNA > 200, percent.mt < 10) %>%
  mutate(patient = str_split(orig.ident, "_", simplify = T)[, 1]) %>%
  mutate(barcode = str_split(cell_id, "_|-", simplify = T)[, 3]) %>%
  mutate(replicate = str_split(cell_id, "-", simplify = T)[, 2])

cell_comp <- cells %>%
  count(patient) %>%
  rename(tomas_nCells = n) %>%
  inner_join(
    geo_cells %>%
      count(patient) %>%
      rename(andrew_nCells = n)
  ) %>%
  mutate(
    added_cells = tomas_nCells - andrew_nCells,
    percent_gain = added_cells / andrew_nCells * 100
  )
print(cell_comp)


df <- cells %>%
  left_join(geo_cells, by = c("patient","barcode"), suffix = c(".tomas", ".andrew"))

df %>% 
  filter(is.na(cell_id.andrew)) %>% 
  count(broad_label)
