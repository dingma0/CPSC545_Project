#!/bin/bash

Rscript /Users/dingma/Desktop/Numbat/scripts/numbat/get_binned_rna.R \
  --rnaCountsFile /Users/dingma/Desktop/Numbat/data/rna/BF15_0736_ATAC_GEX_all_results.rds \
  --outFile /Users/dingma/Desktop/Numbat/data/rna/BF150736_rna_bin.rds \
  --barcodesKeep /Users/dingma/Desktop/Numbat/data/rna/BF150736_rna_mal_barcodes.tsv \
  --geneBinMapCSVFile /Users/dingma/Desktop/Numbat/data/ref/gene2bin_map.csv