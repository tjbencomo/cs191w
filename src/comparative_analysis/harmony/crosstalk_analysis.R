library(readr)
library(dplyr)
library(tidyr)
library(stringr)

pni_cpdb <- read_tsv("data/cpdb/pni-res/pvalues.txt")
ji_cpdb <- read_tsv("data/cpdb/ji-res/pvalues.txt")

info_cols <- colnames(pni_cpdb)[1:11]
pni_df <- pni_cpdb %>% pivot_longer(!all_of(info_cols), names_to = "celltypes", values_to = "pval")
ji_df <- ji_cpdb %>% pivot_longer(!all_of(info_cols), names_to = "celltypes", values_to = "pval")

alpha_thresh <- .01
print(sum(pni_df$pval < alpha_thresh))
print(sum(ji_df$pval < alpha_thresh))

pni_res <- pni_df %>% 
  filter(pval < alpha_thresh) %>%
  mutate(id = str_c(interacting_pair, celltypes, sep=":"))
ji_res <- ji_df %>% 
  filter(pval < alpha_thresh) %>%
  mutate(id = str_c(interacting_pair, celltypes, sep=":"))
