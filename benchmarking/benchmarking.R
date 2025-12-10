# data prep
vec <- c("rna_bin", "rna_gene", "atac_bin", "combined_suf", "combined_sum")
for (mode in vec) {
  segs <- read.table(glue("data/{mode}/segs_consensus_2.tsv"), sep = '\t', header = T) %>%
    select(seg = seg_cons, chr = CHROM, start = seg_start, end = seg_end, length = seg_length, state = cnv_state_post)
  write.csv(segs, glue("results/{mode}/segments.csv"), quote = F, row.names = F)
  
  clones <- read.table(glue("data/{mode}/clone_post_2.tsv"), sep = '\t', header = T) %>%
    select(cell, clone = clone_opt)
  joint <- read.table(glue("data/{mode}/joint_post_2.tsv"), sep = '\t', header = T) %>%
    select(cell, seg, cnv_state_map)
  joined <- full_join(joint, clones, by = 'cell')
  write.csv(joined, glue("results/{mode}/cells.csv"), quote = F, row.names = F)
}


majority_state <- function(x) {
  tab <- table(x)
  if (is.na(tab['neu']) | tab['neu']/sum(tab) < 0.9) {
    names(which.max(tab[c('amp', 'del')]))
  } else {
    'neu'
  }
}

compute_prf <- function(truth, pred, positive_level) {
  truth <- factor(truth, levels = c("neu", "amp", "del"))
  pred  <- factor(pred,  levels = levels(truth))
  
  tp <- sum(truth == positive_level & pred == positive_level, na.rm = TRUE)
  fp <- sum(truth != positive_level & pred == positive_level, na.rm = TRUE)
  fn <- sum(truth == positive_level & pred != positive_level, na.rm = TRUE)
  
  precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
  f1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
               NA, 2 * precision * recall / (precision + recall))
  
  data.frame(class = positive_level,
             tp = tp, fp = fp, fn = fn,
             precision = precision,
             recall = recall,
             f1 = f1)
}

# clean dlp
dlp <- read.table("data/ALL_RESULTS/HMMcopyResults.tsv", header = T) %>% 
  dplyr::filter(quality > 0.75, chr != "chrX", chr != "chrY") %>%
  mutate(cnv_state = case_when(
          state <= 1 ~ 'del',
          state == 2 ~ 'neu',
          state > 2 ~ 'amp'
        ),
        chr = as.integer(str_remove(chr, 'chr'))
  ) %>%
  select(-c('state', 'quality', 'cell_id'))

# dlp bulk summary
dlp_summary <- dlp %>%
  group_by(chr, start, end) %>%
  summarise(
    n_cells = n(),
    frac_amp = mean(cnv_state == "amp"),
    frac_del = mean(cnv_state == "del"),
    frac_neu = mean(cnv_state == "neu"),
    maj_state = majority_state(cnv_state),
    .groups = "drop"
  )

bins_df <- dlp %>%
  distinct(chr, start, end) %>%
  arrange(chr, start)
  
bins_gr <- makeGRangesFromDataFrame(
  bins_df,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = FALSE
)

stats <- data.frame(
  mode = character(),
  accuracy = numeric(),
  cor_amp = numeric(),
  cor_del = numeric(),
  stringsAsFactors = FALSE
)
modes <- c("rna_bin", "rna_gene", "atac_bin", "combined_suf", "combined_sum")
modes <- c("rna_bin")
for (mode in modes) {
  print(glue("{mode}: starting..."))
  # segs to bin numbat
  # clean numbat
  segs <- read.csv(glue("results/{mode}/segments.csv")) %>%
    select(-c('length', 'state'))
  numbat <- read.csv(glue("results/{mode}/cells.csv")) %>%
    select(-clone) %>%
    left_join(segs, by = 'seg') %>%
    mutate(cnv_state_map = case_when(
      cnv_state_map == "bamp" ~ "amp",
      cnv_state_map == "bdel" ~ "del",
      cnv_state_map == "loh" ~ "neu",
      TRUE ~ cnv_state_map
    )) %>%
    select(cell, chr, start, end, cnv_state = cnv_state_map)
  
  seg_gr <- makeGRangesFromDataFrame(
    numbat,
    seqnames.field = "chr",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE,
    na.rm = T
  )
  print(glue("{mode}: creating numbat bins..."))
  ov <- findOverlaps(bins_gr, seg_gr)
  numbat_bin <- as.data.frame(ov) %>%
    transmute(
      chr   = as.character(seqnames(bins_gr[queryHits])),
      start = start(bins_gr[queryHits]),
      end   = end(bins_gr[queryHits]),
      cell = seg_gr$cell[subjectHits],
      cnv_state  = seg_gr$cnv_state[subjectHits]
    ) 
  
  # %>%
  #   group_by(chr, start, end, cell) %>% 
  #   summarise(cnv_state = majority_state(cnv_state), .groups = "drop")
  
  print(glue("{mode}: summarizing..."))
  numbat_bin_summary <- numbat_bin %>%
    group_by(chr, start, end) %>%
    summarise(
      n_cells   = n(),
      frac_amp  = mean(cnv_state == "amp"),
      frac_del  = mean(cnv_state == "del"),
      frac_neu  = mean(cnv_state == "neu"),
      maj_state = majority_state(cnv_state),
      .groups = "drop"
    )
  numbat_bin_summary$chr = as.integer(numbat_bin_summary$chr)
  
  bin_summary <- full_join(
    dlp_summary,
    numbat_bin_summary,
    by = c("chr", "start", "end"),
    suffix = c("_dlp", "_numbat")
  ) %>%
    mutate(
      frac_amp_numbat  = replace_na(frac_amp_numbat, 0),
      frac_del_numbat  = replace_na(frac_del_numbat, 0),
      frac_neu_numbat  = replace_na(frac_neu_numbat, 1),
      maj_state_numbat = replace_na(maj_state_numbat, "neu")
    ) %>%
    select(-n_cells_dlp, -n_cells_numbat)
  
  accuracy <- mean(bin_summary$maj_state_dlp == bin_summary$maj_state_numbat)
  table_states <- table(bin_summary$maj_state_dlp, bin_summary$maj_state_numbat)
  cor_amp <- cor(bin_summary$frac_amp_dlp, bin_summary$frac_amp_numbat, use = "complete.obs")
  cor_del <- cor(bin_summary$frac_del_dlp, bin_summary$frac_del_numbat, use = "complete.obs")
  
  prf_results <- bind_rows(
    compute_prf(bin_summary$maj_state_dlp, bin_summary$maj_state_numbat, "amp"),
    compute_prf(bin_summary$maj_state_dlp, bin_summary$maj_state_numbat, "del"),
    compute_prf(bin_summary$maj_state_dlp, bin_summary$maj_state_numbat, "neu")
  )
  
  write.csv(prf_results, glue("results/{mode}/prf.csv"), quote = F, row.names = F)
  write.csv(bin_summary, glue("results/{mode}/bin_summary.csv"), quote = F, row.names = F)

  stats <- rbind(stats, data.frame(mode = mode, accuracy = accuracy, cor_amp = cor_amp, cor_del = cor_del, stringsAsFactors=FALSE))
  print(glue("{mode}: Done."))
}

write.csv(stats, glue("results/stats.csv"), quote = F, row.names = F)


stats <- read.csv("results/stats.csv")
ggplot(stats, aes(x = mode, y = accuracy)) +
  geom_col(fill = "steelblue", width = 0.5) +
  geom_text(aes(label = round(accuracy, 3)),
            vjust = -0.5, size = 4) +
  ylim(0, 1) +
  labs(
    x = "Mode",
    y = "Accuracy"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


prf <- data.frame(
  class = character(),
  tp = integer(),
  fp = integer(),
  fn = integer(),
  precision = numeric(),
  recall = numeric(),
  f1 = numeric(),
  stringsAsFactors = F
)
modes <- c("rna_bin", "rna_gene", "atac_bin", "combined_suf", "combined_sum")
for (mode in modes) {
  data <- read.csv(glue("results/{mode}/prf.csv")) %>%
    mutate(mode = mode)
  prf <- rbind(prf, data)
}

prf_long <- prf %>%
  pivot_longer(
    cols = c(precision, recall, f1),
    names_to = "metric",
    values_to = "value"
  )

ggplot(prf_long, aes(x = mode, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(value, 2)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.2) +
  facet_wrap(~ class, ncol = 1) +
  scale_fill_manual(values = c(precision = "#1f77b4", recall = "#ff7f0e", f1 = "#2ca02c")) +
  ylim(0, 1) +
  labs(x = "Mode", y = "Metric value", fill = "Metric") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(size = 11, face = "bold")
  ) + scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.2))
  )


# Plot fraction differences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modes <- c("rna_bin", "rna_gene", "atac_bin", "combined_suf", "combined_sum")
plots <- list()
for (mode in modes) {
  bin_summary <- read.csv(glue('results/{mode}/bin_summary.csv'))
  bin_summary <- bin_summary %>%
    mutate(
      mid_mb = (start + end) / 2 / 1e6,
      L1 = abs(frac_amp_dlp - frac_amp_numbat) +
        abs(frac_del_dlp - frac_del_numbat) +
        abs(frac_neu_dlp - frac_neu_numbat)
    )
  
  chr_order <- c(1:22)
  bin_summary2 <- bin_summary %>%
    mutate(chr = factor(chr, levels = chr_order),
           start = as.numeric(start),
           end = as.numeric(end)) %>%
    arrange(chr, start)
  
  chr_offsets <- bin_summary2 %>%
    group_by(chr) %>%
    summarise(chr_len = max(end), .groups = "drop") %>%
    arrange(chr) %>%
    mutate(offset = dplyr::lag(cumsum(chr_len), 0))
  
  bin_summary2 <- bin_summary2 %>%
    left_join(chr_offsets, by = "chr") %>%
    mutate(
      mid_bp = (start + end) / 2,
      genome_pos = mid_bp + offset
    )
  
  chr_centers <- chr_offsets %>%
    mutate(center = offset + chr_len / 2)
  
  plot <- ggplot(bin_summary2, aes(x = genome_pos, y = L1)) +
    geom_line() +
    geom_vline(data = chr_offsets,
               aes(xintercept = offset),
               linetype = "dashed", color = "grey80", linewidth = 0.3) +
    scale_x_continuous(
      breaks = chr_centers$center,
      labels = chr_centers$chr,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Chromosome", y = "L1 distance") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8), panel.spacing = unit(0.5, "lines"))
  plots <- append(plots, plot)
  
  ggsave(glue("results/{mode}/frac_l1.png"), plot, height = 500, width = 3500, units = 'px')
}


combined_plot <- wrap_plots(plots, ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

combined_plot
