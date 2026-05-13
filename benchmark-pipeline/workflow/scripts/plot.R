#!/usr/bin/env Rscript
# Plotting for the chelae benchmark pipeline.
#
# Consumes the merged bench.tsv (and optionally accuracy.tsv) and writes a
# standard bundle of plots into the output dir. Kept deliberately shallow —
# this is "core outputs" only; bespoke exploratory plots belong in notebooks.
#
# Usage:
#   Rscript plot.R <bench_tsv> <out_dir> [accuracy_tsv]

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
bench_tsv <- args[1]
out_dir   <- args[2]
acc_tsv   <- if (length(args) >= 3) args[3] else NA

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bench <- read_tsv(bench_tsv, show_col_types = FALSE) %>%
  mutate(threads = as.integer(threads))

# Aggregate across replicates: median + min/max wall time.
bench_agg <- bench %>%
  group_by(sample, trim_config, tool, threads,
           host_hostname, host_arch, host_cpu_model) %>%
  summarise(
    wall_s_median   = median(wall_s),
    wall_s_min      = min(wall_s),
    wall_s_max      = max(wall_s),
    reads_per_s     = median(reads_per_s),
    bases_per_s     = median(bases_per_s),
    max_rss_kb      = max(max_rss_kb),
    .groups = "drop"
  )

# Throughput vs threads, faceted by (sample × trim_config), colored by tool.
p_throughput <- ggplot(bench_agg,
    aes(x = threads, y = bases_per_s / 1e6, color = tool, group = tool)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = unique(bench_agg$threads)) +
  scale_y_continuous(labels = comma) +
  facet_grid(sample ~ trim_config) +
  labs(title = "Throughput vs thread count",
       x = "threads", y = "bases/sec (millions)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "throughput_vs_threads.pdf"),
       p_throughput, width = 10, height = 7)

# Wall time vs threads, same facets (complement — easier to eyeball).
p_wall <- ggplot(bench_agg,
    aes(x = threads, y = wall_s_median, color = tool, group = tool)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = wall_s_min, ymax = wall_s_max), width = 0.1) +
  scale_x_continuous(breaks = unique(bench_agg$threads)) +
  facet_grid(sample ~ trim_config) +
  labs(title = "Wall time vs thread count",
       x = "threads", y = "seconds (median; bars = min/max)") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "wall_vs_threads.pdf"),
       p_wall, width = 10, height = 7)

# Speedup vs threads — wall_t1 / wall_tN per (sample, trim_config, tool).
# Only meaningful when the run includes threads=1; otherwise emit a
# placeholder so snakemake's declared output is always present.
speedup_path <- file.path(out_dir, "speedup_vs_threads.pdf")
if (1L %in% bench_agg$threads) {
  baseline <- bench_agg %>%
    filter(threads == 1L) %>%
    select(sample, trim_config, tool, wall_t1 = wall_s_median)
  speedup_df <- bench_agg %>%
    inner_join(baseline, by = c("sample", "trim_config", "tool")) %>%
    mutate(speedup = wall_t1 / wall_s_median)
  max_t <- max(speedup_df$threads)
  p_speedup <- ggplot(speedup_df,
      aes(x = threads, y = speedup, color = tool, group = tool)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_line() + geom_point() +
    scale_x_continuous(breaks = unique(speedup_df$threads)) +
    coord_cartesian(xlim = c(1, max_t), ylim = c(1, max_t)) +
    facet_grid(sample ~ trim_config) +
    labs(title = "Parallel speedup (dashed = ideal scaling)",
         x = "threads", y = "speedup vs 1 thread") +
    theme_bw() + theme(legend.position = "bottom")
  ggsave(speedup_path, p_speedup, width = 10, height = 7)
} else {
  placeholder <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = "no speedup plot\n(run did not include threads=1 baseline)") +
    theme_void()
  ggsave(speedup_path, placeholder, width = 6, height = 4)
}

# Max RSS by tool, at default thread count.
p_rss <- bench_agg %>%
  group_by(tool, sample, trim_config) %>%
  summarise(max_rss_kb = max(max_rss_kb), .groups = "drop") %>%
  ggplot(aes(x = tool, y = max_rss_kb / 1024, fill = tool)) +
    geom_col() + coord_flip() +
    facet_grid(sample ~ trim_config) +
    labs(title = "Peak resident memory", x = NULL, y = "max RSS (MB)") +
    theme_bw() + theme(legend.position = "none")
ggsave(file.path(out_dir, "max_rss.pdf"), p_rss, width = 8, height = 6)

# Accuracy heatmap. Always written so snakemake can declare it as output;
# falls back to a "no data" placeholder if the TSV is missing or empty.
acc_path <- file.path(out_dir, "accuracy_heatmap.pdf")
acc_written <- FALSE
if (!is.na(acc_tsv) && file.exists(acc_tsv)) {
  acc <- read_tsv(acc_tsv, show_col_types = FALSE)
  if (nrow(acc) > 0) {
    p_acc <- acc %>%
      filter(trim_config == "adapter_only", mate == "r1", dropped == 0) %>%
      ggplot(aes(x = expected_trim_len, y = observed_trim_len, fill = count)) +
        geom_tile() +
        scale_fill_viridis_c(trans = "log10", na.value = "white") +
        facet_wrap(~ tool) +
        labs(title = "Adapter-trim accuracy (R1, kept reads)",
             x = "expected trim length (bp)", y = "observed trim length (bp)") +
        theme_bw()
    ggsave(acc_path, p_acc, width = 12, height = 8)
    acc_written <- TRUE
  }
}
if (!acc_written) {
  placeholder <- ggplot() +
    annotate("text", x = 0, y = 0, label = "no accuracy data\n(run included no adapter_only trim_configs)") +
    theme_void()
  ggsave(acc_path, placeholder, width = 6, height = 4)
}

cat("Plots written to", out_dir, "\n")
