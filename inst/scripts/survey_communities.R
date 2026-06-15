# Survey plot: loop communities + CTCF + H3K27ac signal tracks -----------
# Author:      JP Flores
# Date:        2026-05-21
# Project:     loopcity
# Description: Addresses Reviewer 2, Comment 2. Produces a plotgardener figure
#              showing loopcity community calls integrated with CTCF and H3K27ac
#              ChIP-seq signal tracks from the same K562 dataset used throughout
#              the paper (Bond et al. 2023, MEGA). Region chr1:110-115 Mb matches
#              Figure 1b/c. Demonstrates how loopcity community calls integrate
#              with other genomic data tracks via plotgardener.
# Input:       data/processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic
#              data/processed/hic/5kbLoops_0.txt (K562 loop calls)
#              data/processed/chip/GSE213908_MERGE_MEGA_ChIP_CTCF_K562_WT_PMA_0_S.bw
#              data/processed/chip/GSE213908_MERGE_MEGA_ChIP_H3K27_K562_WT_PMA_0_S.bw
# Output:      figures/survey_communities.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir  <- "/work/users/j/p/jpflores/projects/loopcity"
data_dir     <- file.path(project_dir, "data")
output_dir   <- file.path(project_dir, "figures")

## Hi-C
hic_file <- file.path(data_dir,
                      "processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic")

## K562 loop calls from Sarah
loops_file <- file.path(data_dir, "processed/hic/5kbLoops_0.txt")

## ChIP-seq bigWigs (0h = undifferentiated K562 control condition)
ctcf_bw <- file.path(data_dir,
                     "processed/chip/GSE213908_MERGE_MEGA_ChIP_CTCF_K562_WT_PMA_0_S.bw")
h3k27ac_bw <- file.path(data_dir,
                        "processed/chip/GSE213908_MERGE_MEGA_ChIP_H3K27_K562_WT_PMA_0_S.bw")

## Communities cache — saved after first run to avoid re-running pipeline
communities_rds <- file.path(data_dir,
                             "processed/loopcity_communities.rds")

## Genomic region — matches Figure 1b/c in the paper
hic_chrom         <- "1"      ## chromosome name in .hic file
bw_chrom          <- "chr1"   ## chromosome name in bigWig files
survey_chromstart <- 110e6
survey_chromend   <- 115e6

## Hi-C display
hic_norm       <- "SCALE"
hic_resolution <- 10e3
hic_zrange     <- c(0, 225)   ## matches Figure 1b/c

## CTCF signal range — capped at 99th percentile to avoid outlier peak
## drowning out the rest of the signal; adjust ctcf_cap as needed
ctcf_cap <- 200

## Page layout (inches)
page_width  <- 9
page_height <- 7.5

## Track heights (inches)
hic_height    <- 2.2
arch_height   <- 0.75
signal_height <- 0.45
gene_height   <- 0.5
gap           <- 0.06

## Colors
ctcf_color      <- "#253494"
h3k27ac_color   <- "#807DBA"
gray_color      <- "#666666"
highlight_color <- "#F0F0F0"   ## very light grey for anchor highlights


# Libraries ---------------------------------------------------------------

library(data.table)
library(GenomeInfoDb)
library(glue)
library(InteractionSet)
library(loopcity)
library(mariner)
library(plotgardener)


# Load data ---------------------------------------------------------------

## Run loopcity pipeline or load cached result.
## Delete communities_rds to force a fresh run.
if (file.exists(communities_rds)) {
  message("Loading cached communities from: ", communities_rds)
  communities <- readRDS(communities_rds)
  message("Seqlevels in cache: ",
          paste(seqlevels(communities), collapse = ", "))
} else {
  merged <- data.table::fread(loops_file) |>
    mariner::assignToBins(binSize = 10e3) |>
    mergeAnchors(pixelOverlap = 1)
  GenomeInfoDb::seqlevelsStyle(merged) <- "ENSEMBL"
  connected   <- connectLoopAnchors(merged, overlapDist = 1e6)
  scored      <- scoreInteractions(connected,
                                   hicFile    = hic_file,
                                   resolution = hic_resolution,
                                   norm       = hic_norm)
  communities <- assignCommunities(interactions(scored),
                                   leidenResolution = 0.5)
  saveRDS(communities, communities_rds)
}


# Wrangle -----------------------------------------------------------------

## Filter to loops where both anchors fall within the survey region,
## have a valid community assignment, and score > 0
survey_region <- GenomicRanges::GRanges(
  seqnames = hic_chrom,
  ranges   = IRanges::IRanges(start = survey_chromstart,
                              end   = survey_chromend)
)

a1_in <- GenomicRanges::countOverlaps(
  InteractionSet::anchors(communities, "first"),  survey_region) > 0
a2_in <- GenomicRanges::countOverlaps(
  InteractionSet::anchors(communities, "second"), survey_region) > 0

loops_region <- communities[a1_in & a2_in]
loops_region <- loops_region[
  which(loops_region$score > 0 &
          lengths(loops_region$loopCommunity) > 0)
]

## loopcityColors is an internal function — define it here directly
loopcity_colors <- function(n) {
  base_cols <- c("chartreuse3", "deepskyblue3", "darkorange",
                 "darkorchid2", "deeppink3")
  rep(base_cols, length.out = n)
}

## Extract first community label per loop; NA if unassigned
comm_labels <- vapply(
  loops_region$loopCommunity,
  function(x) {
    if (length(x) == 0L) NA_integer_
    else as.integer(x[[1L]])
  },
  integer(1L)
)

message(glue("Community labels: {sum(!is.na(comm_labels))} assigned, ",
             "{sum(is.na(comm_labels))} unassigned"))
message(glue("Unique communities: ",
             "{paste(sort(unique(comm_labels)), collapse = ', ')}"))

## Re-index to 1-based for color palette
comm_labels_reindexed <- as.integer(factor(comm_labels))
n_comms               <- max(comm_labels_reindexed, na.rm = TRUE)
loops_region$color    <- loopcity_colors(n_comms)[comm_labels_reindexed]
loops_region$color[is.na(comm_labels)] <- "grey70"

## Coerce DelayedMatrix to plain numeric
loops_region$score_lg <- as.numeric(log2(loops_region$score))

## Get unique merged anchors in the survey region for highlights
## Both anchor1 and anchor2 positions, deduplicated
anchor1_gr <- InteractionSet::anchors(loops_region, "first")
anchor2_gr <- InteractionSet::anchors(loops_region, "second")
all_anchors <- sort(unique(c(anchor1_gr, anchor2_gr)))


# Diagnostics -------------------------------------------------------------

message(glue("loops_region has {length(loops_region)} loops"))
message(glue("score_lg range: ",
             "{round(min(loops_region$score_lg, na.rm=TRUE), 2)} ",
             "to {round(max(loops_region$score_lg, na.rm=TRUE), 2)}"))
message(glue("Unique anchors to highlight: {length(all_anchors)}"))


# Visualization -----------------------------------------------------------

## Shared pgParams — all tracks lock to the same coordinates
params <- pgParams(
  assembly   = "hg38",
  chrom      = bw_chrom,
  chromstart = survey_chromstart,
  chromend   = survey_chromend,
  resolution = hic_resolution,
  zrange     = hic_zrange,
  norm       = hic_norm
)

## Signal ranges
## CTCF: capped at ctcf_cap to prevent outlier peak from drowning the track
ctcf_range    <- c(0, ctcf_cap)
h3k27ac_range <- calcSignalRange(
  h3k27ac_bw,
  chrom      = bw_chrom,
  chromstart = survey_chromstart,
  chromend   = survey_chromend,
  assembly   = "hg38",
  negData    = FALSE
)

## y positions computed top-to-bottom
x     <- 0.5
width <- page_width - 1.0

y_hic     <- 0.4
y_arch    <- y_hic     + hic_height    + gap
y_ctcf    <- y_arch    + arch_height   + gap * 0.25
y_h3k27ac <- y_ctcf    + signal_height + gap * 3
y_genes   <- y_h3k27ac + signal_height + gap * 0.5
y_label   <- y_genes   + gene_height   + 0.03

## Highlights span from bottom of Hi-C to top of gene track
highlight_y      <- y_hic + hic_height
highlight_height <- (y_h3k27ac + signal_height) - highlight_y


# Save outputs ------------------------------------------------------------

pdf(file   = file.path(output_dir, "survey_communities.pdf"),
    width  = page_width,
    height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

## ── Hi-C contact map ──────────────────────────────────────────────────────
hic_plot <- plotHicRectangle(
  data   = hic_file,
  params = params,
  chrom  = hic_chrom,
  zrange = hic_zrange,
  x      = x,
  y      = y_hic,
  width  = width,
  height = hic_height
)

annoHeatmapLegend(
  hic_plot,
  x      = x - 0.25,
  y      = y_hic,
  width  = 0.07,
  height = hic_height / 3,
  just   = c("left", "top"),
  default.units = "inches"
)

## ── Anchor highlights — light grey boxes spanning all tracks ──────────────
for (i in seq_along(all_anchors)) {
  annoHighlight(
    plot       = hic_plot,
    chrom      = hic_chrom,
    chromstart = start(all_anchors[i]),
    chromend   = end(all_anchors[i]),
    y          = highlight_y,
    height     = highlight_height,
    just       = c("left", "top"),
    fill       = highlight_color,
    alpha      = 0.5,
    default.units = "inches"
  )
}

## ── Community-coloured loop arches ────────────────────────────────────────
if (length(loops_region) > 0) {
  score_range <- c(min(loops_region$score_lg, na.rm = TRUE),
                   max(loops_region$score_lg, na.rm = TRUE))
  plotPairsArches(
    loops_region,
    params     = params,
    chrom      = hic_chrom,
    flip       = TRUE,
    clip       = TRUE,
    archHeight = "score_lg",
    range      = score_range,
    x          = x,
    y          = y_arch,
    width      = width,
    height     = arch_height,
    fill       = loops_region$color,
    linecolor  = loops_region$color
  )
}

## ── CTCF ChIP-seq signal ──────────────────────────────────────────────────
plotSignal(
  data      = ctcf_bw,
  params    = params,
  chrom     = bw_chrom,
  x         = x,
  y         = y_ctcf,
  width     = width,
  height    = signal_height,
  fill      = ctcf_color,
  linecolor = ctcf_color,
  range     = ctcf_range,
  scale     = FALSE
)
plotText(
  label     = "CTCF",
  x         = x,
  y         = y_ctcf - 0.05,
  fontsize  = 7,
  fontface  = "bold",
  fontcolor = ctcf_color,
  just      = c("left", "bottom")
)
plotText(
  label     = glue("0-{ctcf_cap}"),
  x         = x + width,
  y         = y_ctcf + 0.03,
  fontsize  = 5,
  fontcolor = gray_color,
  just      = c("right", "top")
)

## ── H3K27ac ChIP-seq signal ───────────────────────────────────────────────
plotSignal(
  data      = h3k27ac_bw,
  params    = params,
  chrom     = bw_chrom,
  x         = x,
  y         = y_h3k27ac,
  width     = width,
  height    = signal_height,
  fill      = h3k27ac_color,
  linecolor = h3k27ac_color,
  range     = h3k27ac_range,
  scale     = FALSE
)
plotText(
  label     = "H3K27ac",
  x         = x,
  y         = y_h3k27ac - 0.05,
  fontsize  = 7,
  fontface  = "bold",
  fontcolor = h3k27ac_color,
  just      = c("left", "bottom")
)
plotText(
  label     = glue("0-{ceiling(h3k27ac_range[2])}"),
  x         = x + width,
  y         = y_h3k27ac + 0.03,
  fontsize  = 5,
  fontcolor = gray_color,
  just      = c("right", "top")
)

## ── Gene models ───────────────────────────────────────────────────────────
plotGenes(
  params   = params,
  x        = x,
  y        = y_genes,
  width    = width,
  height   = gene_height,
  fontsize = 6
)

## ── Genome coordinate label ───────────────────────────────────────────────
annoGenomeLabel(
  plot     = hic_plot,
  x        = x,
  y        = y_label,
  scale    = "bp",
  fontsize = 6
)

dev.off()


# Session info ------------------------------------------------------------

sessionInfo()