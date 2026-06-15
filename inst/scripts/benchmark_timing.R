# Benchmark loopcity runtime — full pipeline, whole genome, 100 iters ----
# Author:      JP Flores
# Date:        2026-05-27
# Project:     loopcity
# Description: Times each step of the loopcity pipeline independently over
#              100 iterations using bench::mark(). Each step's input is
#              pre-computed once so only the target function is timed.
#              Produces a horizontal bar plot (median ± IQR) with function
#              names on the y-axis and time in seconds on the x-axis.
#              Uses K562 chr21 loops and Hi-C file (Bond et al. 2023),
#              10 kb resolution.
# Input:       data/processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic
#              data/processed/hic/5kbLoops_0.txt (K562 loop calls)
# Output:      man/figures/benchmark_timing.pdf
#              data/processed/benchmark_timing.rds
#              prints manuscript sentence to console
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir    <- "/work/users/j/p/jpflores/projects/loopcity"
data_dir       <- file.path(project_dir, "data")
output_data    <- file.path(project_dir, "data/processed")

hic_file       <- file.path(data_dir,
                            "processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic")
loops_file     <- file.path(data_dir, "processed/hic/5kbLoops_0.txt")

n_iterations   <- 100
hic_resolution <- 10e3
hic_norm       <- "NONE"
pixel_overlap  <- 1
overlap_dist   <- 1e6
leiden_res     <- 0.5


# Libraries ---------------------------------------------------------------

library(bench)
library(data.table)
library(GenomeInfoDb)
library(ggplot2)
library(glue)
library(InteractionSet)
library(loopcity)
library(mariner)


# Load data ---------------------------------------------------------------

all_loops <- data.table::fread(loops_file) |>
  mariner::assignToBins(binSize = 10e3)
GenomeInfoDb::seqlevelsStyle(all_loops) <- "ENSEMBL"

## Subset to chr21 — smallest autosome, representative for benchmarking
loops <- all_loops[GenomeInfoDb::seqnames(
  InteractionSet::anchors(all_loops, "first")) == "21"]

message(glue("chr21 loops: {length(loops)}"))


# Pre-compute inputs for each step ----------------------------------------

## Each step's input is computed once so bench::mark() only times
## the target function, not its upstream dependencies.
merged    <- mergeAnchors(loops, pixelOverlap = pixel_overlap)
connected <- connectLoopAnchors(merged, overlapDist = overlap_dist)
scored    <- scoreInteractions(connected,
                               hicFile    = hic_file,
                               resolution = hic_resolution,
                               norm       = hic_norm)


# Benchmark ---------------------------------------------------------------

message(glue("Running {n_iterations} iterations per step..."))

bm_merge <- bench::mark(
  mergeAnchors(loops, pixelOverlap = pixel_overlap),
  iterations = n_iterations, check = FALSE, memory = FALSE
)

bm_connect <- bench::mark(
  connectLoopAnchors(merged, overlapDist = overlap_dist),
  iterations = n_iterations, check = FALSE, memory = FALSE
)

bm_score <- bench::mark(
  scoreInteractions(connected,
                    hicFile    = hic_file,
                    resolution = hic_resolution,
                    norm       = hic_norm),
  iterations = n_iterations, check = FALSE, memory = FALSE
)

bm_community <- bench::mark(
  assignCommunities(interactions(scored),
                    leidenResolution = leiden_res),
  iterations = n_iterations, check = FALSE, memory = FALSE
)


# Wrangle -----------------------------------------------------------------

extract_times <- function(bm, step_label) {
  times <- as.numeric(bm$time[[1]])
  data.frame(
    step     = step_label,
    median_s = median(times),
    q25_s    = quantile(times, 0.25),
    q75_s    = quantile(times, 0.75)
  )
}

step_times <- rbind(
  extract_times(bm_merge,     "mergeAnchors"),
  extract_times(bm_connect,   "connectLoopAnchors"),
  extract_times(bm_score,     "scoreInteractions"),
  extract_times(bm_community, "assignCommunities")
)

## Synthetic Total row — sum of per-step medians and IQR bounds
total_row <- data.frame(
  step     = "Total",
  median_s = sum(step_times$median_s),
  q25_s    = sum(step_times$q25_s),
  q75_s    = sum(step_times$q75_s)
)

timing <- rbind(step_times, total_row)
timing$step <- factor(timing$step,
                      levels = c("Total", "assignCommunities",
                                 "scoreInteractions", "connectLoopAnchors",
                                 "mergeAnchors"))

message("\nPipeline timing summary (median seconds):")
print(timing)

## Manuscript sentence
total_min  <- round(timing$median_s[timing$step == "Total"] / 60, 2)
slow_step  <- timing$step[which.max(
  timing$median_s[timing$step != "Total"])]
slow_min   <- round(timing$median_s[timing$step == slow_step] / 60, 2)

message(glue(
  "\nManuscript sentence:",
  "\nThe full loopcity pipeline completed in {total_min} minutes on the ",
  "complete K562 loop call set (whole genome, 10 kb resolution), with ",
  "{slow_step} representing the most time-intensive step ",
  "({slow_min} minutes)."
))


# Visualize ---------------------------------------------------------------

step_labels <- c(
  mergeAnchors       = "mergeAnchors()",
  connectLoopAnchors = "connectLoopAnchors()",
  scoreInteractions  = "scoreInteractions()",
  assignCommunities  = "assignCommunities()",
  Total              = "Total"
)

bar_fill <- ifelse(timing$step == "Total", "#2166AC", "#4393C3")

p_timing <- ggplot(timing, aes(x = median_s, y = step)) +
  geom_col(fill = bar_fill, width = 0.65) +
  geom_errorbar(aes(xmin = q25_s, xmax = q75_s),
                width = 0.25, linewidth = 0.4, color = "grey30") +
  geom_text(
    data = ~subset(.x, !step %in% c("scoreInteractions", "Total")),
    aes(label = round(median_s, 1)),
    hjust = -0.2, size = 3, family = "mono"
  ) +
  geom_text(
    data = ~subset(.x, step %in% c("scoreInteractions", "Total")),
    aes(label = round(median_s, 1)),
    hjust = -0.2, size = 3, family = "mono", nudge_y = 0.35
  ) +
  scale_y_discrete(labels = step_labels) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Time (seconds)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y  = element_text(size = 9, family = "mono"),
    axis.text.x  = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    axis.line    = element_line(linewidth = 0.3),
    axis.ticks   = element_line(linewidth = 0.3)
  )


# Whole-genome single run for manuscript sentence -------------------------

message("Running whole-genome pipeline once for manuscript sentence...")

t_total <- system.time({
  merged_all      <- mergeAnchors(all_loops, pixelOverlap = pixel_overlap)
  connected_all   <- connectLoopAnchors(merged_all, overlapDist = overlap_dist)
  scored_all      <- scoreInteractions(connected_all,
                                       hicFile    = hic_file,
                                       resolution = hic_resolution,
                                       norm       = hic_norm)
  communities_all <- assignCommunities(interactions(scored_all),
                                       leidenResolution = leiden_res)
})

total_min_wg <- round(t_total[["elapsed"]] / 60, 2)

message(glue(
  "
Manuscript sentence:",
  "
The full loopcity pipeline completed in {total_min_wg} minutes on the ",
  "complete K562 loop call set (whole genome, 10 kb resolution)."
))

# Save outputs ------------------------------------------------------------

saveRDS(timing, file = file.path(output_data, "benchmark_timing.rds"))

ggsave(
  filename = file.path(project_dir, "man/figures/benchmark_timing.pdf"),
  plot     = p_timing,
  width    = 5,
  height   = 4
)


# Session info ------------------------------------------------------------

sessionInfo()