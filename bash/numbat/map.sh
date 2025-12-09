#!/bin/bash

Rscript /Users/dingma/Desktop/Numbat/scripts/numbat/get_gene_binned_intersections.R \
  --numbatGTFname hg38 \
  --binGR /Users/dingma/Desktop/Numbat/data/ref/var220kb.rds \
  --outfile /Users/dingma/Desktop/Numbat/data/ref/gene2bin_map.csv
  