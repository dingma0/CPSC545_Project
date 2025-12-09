source("scripts/numbat/input_prep.R")

saveRDS(binCnt_union(
  c("data/rna/BF150736_rna_bin.rds",
    "data/atac/BF150736_atac_bin.rds"),
  seed=123,maxCB=10000,
  add_suffix = FALSE, 
  runMultiomeAsSummed = TRUE
), "data/both/sum/BF150736_comb_sum_bincnt.rds")

combine_allele_counts(
  alleleCounts_RNA = "data/rna/BF150736_rna_allele_counts.tsv.gz",
  alleleCounts_ATAC = "data/atac/BF150736_atac_allele_counts.tsv.gz",
  output_file = "data/both/sum/BF150736_comb_sum_allele_counts.tsv.gz",
  addBarcodeSuff = FALSE
)

# ref
saveRDS(binCnt_union(
  c("data/rna/BF150736_norm_rna_bin.rds",
    "data/atac/BF150736_norm_atac_bin.rds"),
  seed=123,maxCB=10000,
  add_suffix = FALSE, 
  runMultiomeAsSummed = TRUE
), "data/both/sum/BF150736_norm_comb_sum_bincnt.rds")
