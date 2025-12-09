#!/bin/bash

Rscript /mnt/scripts/run_numbat_multiome.R  \
            --countmat /mnt/data/rna_bin/BF150736_rna_bin.rds \
            --alleledf /mnt/data/rna_bin/BF150736_rna_allele_counts.tsv.gz \
            --out_dir /mnt/results/BF150736_rna_bin \
            --ref /mnt/data/rna_bin/BF150736_norm_rna_bin.rds \
            --gtf /mnt/data/var220kb.rds \
            --parL /mnt/data/par_numbatm.rds
