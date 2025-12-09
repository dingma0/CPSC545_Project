source("scripts/numbat/input_prep.R")

saveRDS(binCnt(
  c("data/rna/BF150736_rna_bin.rds",
    "data/atac/BF150736_atac_bin.rds"),
  seed=123,maxCB=10000,
  add_suffix = TRUE, 
  suffix_tag = c('_RNA', '_ATAC')
), "data/both/BF150736_comb_bincnt.rds")

combine_allele_counts(
  alleleCounts_RNA = "data/rna/BF150736_rna_allele_counts.tsv.gz",
  alleleCounts_ATAC = "data/atac/BF150736_atac_allele_counts.tsv.gz",
  output_file = "data/both/BF150736_comb_allele_counts.tsv.gz",
  addBarcodeSuff = TRUE
)

ref_rna <- readRDS("data/rna/BF150736_norm_rna_bin.rds")
ref_atac <- readRDS("data/atac/BF150736_norm_atac_bin.rds")
shared <- intersect(rownames(ref_rna), rownames(ref_atac))
ref_comb <- cbind(ref_rna[shared, ], ref_atac[shared, ])
saveRDS(ref_comb, "data/both/BF150736_norm_comb_bincnt.rds")