#!/usr/bin/env Rscript
# plot_lengths_R.R
# Average reinflated length distribution (15-90nt), one line per group,
# from *_ge5_detected.hist files (weighted length<TAB>count, no header).
#
# Ported from the standalone plot_lengths_R_ggplot.R used to build this
# project's real deliverables -- generalised from that script's R*->RJ /
# T*->SLT filename-prefix guess to a real --samplesheet-driven group lookup
# (any number of groups, any labels), same convention as every other script
# in this pipeline. Percentages are computed over each sample's full hist
# range (not windowed before averaging); the plot itself crops the x-axis
# display to 15-90nt.

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(magrittr)
  library(reshape2)
  library(svglite)
  library(data.table)
})

# Simple base-R positional/flag parsing (no argparse dependency):
#   plot_lengths_R.R <hist_dir> --samplesheet <path> --outdir <path>
raw_args <- commandArgs(trailingOnly = TRUE)
flags <- list(); positional <- c()
i <- 1
while (i <= length(raw_args)) {
  a <- raw_args[i]
  if (a %in% c("--samplesheet", "--outdir")) {
    flags[[sub("^--", "", a)]] <- raw_args[i + 1]
    i <- i + 2
  } else {
    positional <- c(positional, a)
    i <- i + 1
  }
}
if (length(positional) < 1) stop("Usage: plot_lengths_R.R <hist_dir> --samplesheet <path> --outdir <path>")
if (is.null(flags$samplesheet)) stop("Missing required argument: --samplesheet")
if (is.null(flags$outdir))      stop("Missing required argument: --outdir")
args <- list(hist_dir = positional[1], samplesheet = flags$samplesheet, outdir = flags$outdir)

load_samplesheet <- function(path) {
  sheet <- read.csv(path, stringsAsFactors = FALSE)
  setNames(as.character(sheet$group), as.character(sheet$sample))
}
sample_to_group <- load_samplesheet(args$samplesheet)

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

hist_files <- list.files(args$hist_dir, pattern = "_ge5_detected\\.hist$", full.names = TRUE)
if (length(hist_files) == 0) {
  stop(sprintf("No *_ge5_detected.hist files found in %s", args$hist_dir))
}

sample_name <- function(path) sub("_ge5_detected\\.hist$", "", basename(path))
infer_group <- function(name) {
  g <- sample_to_group[[name]]
  if (is.null(g)) NA_character_ else g
}

names_all  <- vapply(hist_files, sample_name, character(1))
groups_all <- vapply(names_all, infer_group, character(1))
groups     <- sort(unique(groups_all[!is.na(groups_all)]))

if (length(groups) == 0) {
  stop("No hist files matched any sample in the samplesheet")
}

load_group <- function(files, mean_col) {
  samples <- lapply(files, function(f) as.data.frame(data.table::fread(f)))
  names(samples) <- sample_name(files)

  samples <- lapply(samples, function(x) {
    x[, 2] <- x[, 2] / sum(x[, 2])
    x
  })

  samples <- imap(samples, function(df, name) {
    colnames(df)[1] <- "Length"
    colnames(df)[2] <- name
    df
  })

  merged <- samples %>% purrr::reduce(full_join, by = "Length")
  merged[is.na(merged)] <- 0
  merged[[mean_col]] <- rowMeans(merged[, -1, drop = FALSE])
  merged
}

# Load one averaged curve per group, keeping only Length + that group's mean col.
group_curves <- list()
for (g in groups) {
  files <- hist_files[groups_all == g]
  mean_col <- paste0(g, "_mean")
  merged <- load_group(files, mean_col)
  group_curves[[g]] <- merged[, c("Length", mean_col)]
}

all_together <- Reduce(function(a, b) dplyr::full_join(a, b, by = "Length"), group_curves)
all_together[is.na(all_together)] <- 0

# Zero-anchor at the window edges so every curve visually starts/ends at 0%
# instead of floating at whatever value it has at 15/90nt.
mean_cols <- paste0(groups, "_mean")
zero_row <- function(len) {
  row <- as.list(setNames(rep(0, length(mean_cols)), mean_cols))
  row$Length <- len
  as.data.frame(row)[, c("Length", mean_cols)]
}
all_together <- dplyr::bind_rows(zero_row(15), all_together, zero_row(90))

all_together <- all_together %>% reshape2::melt(id.vars = "Length")

# Cycled through in order for however many groups are present -- same
# palette convention as plot_lengths.py/plot_report_summary.py.
color_pool <- c("#00a4ff", "#ff0000", "#2ca02c", "#9467bd",
                 "#ff7f0e", "#17becf", "#e377c2", "#8c564b")
group_colors <- setNames(color_pool[seq_along(groups)], mean_cols)
group_labels <- setNames(groups, mean_cols)

p <- ggplot(all_together, aes(x = Length, y = value, colour = variable)) +
  geom_line(linewidth = 0.3) +
  theme_bw(base_size = 5, base_family = "Nimbus Sans") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(15, 90),
                      breaks = c(15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)) +
  scale_colour_manual(values = group_colors, labels = group_labels) +
  labs(x = "Sequence length (nt)", y = "Percentage of reads (%)", colour = NULL) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(6, "pt"),
    legend.text = element_text(size = 5),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0, "pt"),
    axis.text = element_text(size = 5),
    panel.grid.major.x = element_line(colour = "grey65", linetype = "dotted"),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor.y = element_line(colour = "grey90")
  )

today <- format(Sys.Date(), "%d_%m_%Y")
outfile_svg <- file.path(args$outdir, paste0(today, "_avg_reinflated_length_distribution_R.svg"))
ggsave(outfile_svg, p, width = 168, height = 98, units = "px", dpi = 72, device = svglite::svglite)
cat("Saved:", outfile_svg, "\n")
