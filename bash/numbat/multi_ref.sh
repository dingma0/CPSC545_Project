#!/bin/bash

Rscript /Users/dingma/Desktop/Numbat/scripts/numbat/get_binned_rna.R \
  --rnaCountsFile /Users/dingma/Desktop/Numbat/data/rna/BF15_0736_ATAC_GEX_all_results.rds \
  --outFile /Users/dingma/Desktop/Numbat/data/rna/BF150736_norm_rna_bin.rds \
  --barcodesKeep /Users/dingma/Desktop/Numbat/data/normal_cells/BF150736_norm_barcodes.tsv \
  --geneBinMapCSVFile /Users/dingma/Desktop/Numbat/data/ref/gene2bin_map.csv \
  --generateAggRef
  
Rscript /Users/dingma/Desktop/Numbat/scripts/numbat/get_binned_atac.R \
  --CB /Users/dingma/Desktop/Numbat/data/normal_cells/BF150736_norm_barcodes.tsv \
  --frag /Users/dingma/Desktop/Numbat/data/atac/atac_fragments.tsv.gz \
  --binGR /Users/dingma/Desktop/Numbat/data/ref/var220kb.rds \
  --outFile /Users/dingma/Desktop/Numbat/data/atac/BF150736_norm_atac_bin.rds \
  --generateAggRef