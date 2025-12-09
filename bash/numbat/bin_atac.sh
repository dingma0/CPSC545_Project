#!/bin/bash

Rscript /Users/dingma/Desktop/Numbat/scripts/numbat/get_binned_atac.R \
  --CB /Users/dingma/Desktop/Numbat/data/atac/BF150736_atac_mal_barcodes.tsv \
  --frag /Users/dingma/Desktop/Numbat/data/atac/atac_fragments.tsv.gz \
  --binGR /Users/dingma/Desktop/Numbat/data/ref/var220kb.rds \
  --outFile /Users/dingma/Desktop/Numbat/data/atac/BF150736_atac_bin.rds