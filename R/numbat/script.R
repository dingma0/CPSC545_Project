norm <- read.csv("data/normal_cells/atac_7_samples_normal_cell_ids.csv")

norm <- dplyr::filter(norm, sample == 'BF15_0736') %>%
  select(cell, group = seurat_clusters) %>%
  mutate(cell = str_remove(cell, "BF15_0736_"))

# write.table(norm, "data/normal_cells/BF150736_norm_barcodes.tsv", sep = '\t', quote = F, row.names = F)

atac_all <- read.csv("data/filtered_feature_bc_matrix/barcodes.tsv.gz", header = F, col.names = "cell")
atac_mal <- setdiff(atac_all$cell, norm$cell)
# write.csv(atac_mal, "results/BF150736_atac_mal_barcodes.tsv", quote = F, row.names = F)

rna_all <- read.csv("data/BF150736_rna_barcodes.tsv", header = F, col.names = 'cell')
rna_mal <- setdiff(rna_all$cell, norm$cell)
# write.csv(rna_mal, "results/BF150736_rna_mal_barcodes.tsv", quote = F, row.names = F)

par_numbatm <- list(
  genome = "hg38",
  t = 1e-5,
  ncores = 32,
  plot = TRUE
)
# saveRDS(par_numbatm, "data/both/par_numbatm.rds")


source("scripts/numbat/input_prep.R")

saveRDS(binCnt_union(
  c("data/rna/BF150736_rna_bin.rds",
    "data/atac/BF150736_atac_bin.rds"),
  seed=123,maxCB=10000,
  add_suffix = TRUE, 
  suffix_tag = c('_RNA', '_ATAC')
), "data/temp/BF150736_comb_union_suf_bincnt.rds")


obj <- readRDS("data/temp/BF150736_comb_union_suf_bincnt.rds")
obj0 <- readRDS("data/both/suf/BF150736_comb_bincnt.rds")
diff <- obj[setdiff(rownames(obj), rownames(obj0)),]



