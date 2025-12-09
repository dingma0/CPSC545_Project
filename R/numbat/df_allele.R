label = 'BF150736'
samples = 'BF150736_rna'
outdir = '/Users/dingma/Desktop/Numbat/data/temp/out'
n_samples = length(samples)
gmap = 'data/temp/genetic_map_hg38_withX.txt.gz'
snpvcf = 'data/temp/BF150736_phased_snps.vcf.gz'
genome = 'hg38'

## VCF creation
cat('Creating VCFs\n')

# read in the pileup VCF
vcfs = lapply(samples, function(sample) {
    vcf_file = glue('data/temp/cellSNP.base.vcf')
    if (file.exists(vcf_file)) {
        if (file.size(vcf_file) != 0) {
            vcf = vcfR::read.vcfR(vcf_file, verbose = F)
            if (nrow(vcf@fix) == 0) {
                stop(glue('Pileup VCF for sample {sample} has 0 variants'))
            }
            return(vcf)
        } else {
            stop('Pileup VCF is empty')
        }
    } else {
        stop('Pileup VCF not found')
    }
})

# Remove chr prefix if present
vcfs = lapply(vcfs, function(vcf){
    vcf@fix[,1] <- gsub("chr", "", vcf@fix[,1])
    return(vcf)
})

## Generate allele count dataframe
cat('Generating allele count dataframes\n')

if (genome == 'hg19') {
    gtf = gtf_hg19
} else {
    gtf = gtf_hg38
}

genetic_map = fread(gmap) %>% 
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

    
# read in phased VCF
vcf_phased <- fread(snpvcf, skip = '#CHROM') %>%
    dplyr::rename(CHROM = `#CHROM`) %>%
    dplyr::mutate(CHROM = str_remove(CHROM, 'chr')) %>%
    mutate(CHROM = factor(CHROM, unique(CHROM)))

last_col <- names(vcf_phased)[ncol(vcf_phased)]
vcf_phased <- mutate(vcf_phased, !!label := stringr::str_remove(get(last_col), ":.*"))

pu_dir = 'data/temp'

# pileup VCF
vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM') %>% 
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, 'chr'))

# count matrices
AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))

cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# pileup VCF
vcf_pu = vcf_pu %>%
  mutate(INFO = str_remove_all(INFO, '[:alpha:]|=')) %>%
  tidyr::separate(col = 'INFO', into = c('AD', 'DP', 'OTH'), sep = ';') %>%
  mutate_at(c('AD', 'DP', 'OTH'), as.integer) %>%
  mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_'))

# pileup count matrices
DP = as.data.frame(Matrix::summary(DP)) %>%
  mutate(
    cell = cell_barcodes[j],
    snp_id = vcf_pu$snp_id[i]
  ) %>%
  select(-i, -j) %>%
  rename(DP = x) %>%
  select(cell, snp_id, DP)

AD = as.data.frame(Matrix::summary(AD)) %>%
  mutate(
    cell = cell_barcodes[j],
    snp_id = vcf_pu$snp_id[i]
  ) %>%
  select(-i, -j) %>%
  rename(AD = x) %>%
  select(cell, snp_id, AD)

df = DP %>% left_join(AD, by = c("cell", "snp_id")) %>%
  mutate(AD = ifelse(is.na(AD), 0, AD))

df = df %>% left_join(
  vcf_pu %>% rename(AD_all = AD, DP_all = DP, OTH_all = OTH),
  by = 'snp_id')

df = df %>% mutate(
  AR = AD/DP,
  AR_all = AD_all/DP_all
)

df = df %>% dplyr::filter(DP_all > 1 & OTH_all == 0)

# vcf has duplicated records sometimes
df = df %>% distinct()

df = df %>% mutate(
  snp_index = as.integer(factor(snp_id, unique(snp_id))),
  cell_index = as.integer(factor(cell, sample(unique(cell))))
)

# phased VCF
vcf_phased = vcf_phased %>% 
  mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_')) %>%
  mutate(GT = get(label))

# annotate SNP by gene
overlap_transcript = GenomicRanges::findOverlaps(
  vcf_phased %>% {GenomicRanges::GRanges(
    seqnames = .$CHROM,
    IRanges::IRanges(start = .$POS,
                     end = .$POS)
  )},
  gtf %>% {GenomicRanges::GRanges(
    seqnames = .$CHROM,
    IRanges::IRanges(start = .$gene_start,
                     end = .$gene_end)
  )}
) %>%
  as.data.frame() %>%
  setNames(c('snp_index', 'gene_index')) %>%
  left_join(
    vcf_phased %>% mutate(snp_index = 1:n()) %>%
      select(snp_index, snp_id),
    by = c('snp_index')
  ) %>%
  left_join(
    gtf %>% mutate(gene_index = 1:n()),
    by = c('gene_index')
  ) %>%
  arrange(snp_index, gene) %>%
  distinct(snp_index, `.keep_all` = TRUE)

vcf_phased = vcf_phased %>%
  left_join(
    overlap_transcript %>% select(snp_id, gene, gene_start, gene_end),
    by = c('snp_id')
  )

# annotate SNPs by genetic map
marker_map = GenomicRanges::findOverlaps(
  vcf_phased %>% {GenomicRanges::GRanges(
    seqnames = .$CHROM,
    IRanges::IRanges(start = .$POS,
                     end = .$POS)
  )},
  genetic_map %>% {GenomicRanges::GRanges(
    seqnames = .$CHROM,
    IRanges::IRanges(start = .$start,
                     end = .$end)
  )}
) %>%
  as.data.frame() %>%
  setNames(c('marker_index', 'map_index')) %>%
  left_join(
    vcf_phased %>% mutate(marker_index = 1:n()) %>%
      select(marker_index, snp_id),
    by = c('marker_index')
  ) %>%
  left_join(
    genetic_map %>% mutate(map_index = 1:n()),
    by = c('map_index')
  ) %>%
  arrange(marker_index, -start) %>%
  distinct(marker_index, `.keep_all` = TRUE) %>%
  select(snp_id, cM)

vcf_phased = vcf_phased %>% 
  left_join(marker_map, by = 'snp_id')

# add annotation to cell counts
df = df %>% left_join(vcf_phased %>% select(snp_id, gene, GT, cM), by = 'snp_id')
df = df %>% mutate(CHROM = factor(CHROM, unique(CHROM)))
df = df %>% select(cell, snp_id, CHROM, POS, cM, REF, ALT, AD, DP, GT, gene)

# only keep hets
df = df %>% filter(GT %in% c('1|0', '0|1'))

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

fwrite(df, glue('{outdir}/{samples}_allele_counts.tsv.gz'), sep = '\t')

cat('All done!\n')