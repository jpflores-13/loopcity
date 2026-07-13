# TAD null-model sparsity comparison ---------------------------------------
# Author:      JP Flores
# Date:        2026-07-13
# Project:     loopcity
# Description: Addresses reviewer comment on comparing loopcity's community
#              network against a naive null model: assume all merged loop
#              anchors within a single TAD form a fully connected graph.
#              Quantifies how much sparser the network becomes after the
#              loopcity community-detection step, demonstrating specificity
#              of called loops versus mere spatial (TAD) proximity.
#
#              TADs come from Arrowhead calls on the same K562 MEGA dataset
#              used throughout the paper (Bond et al. 2023), 0h condition,
#              matching the .hic file used in survey_communities.R and
#              benchmark_timing.R. See:
#              github.com/mbond0718/MEGA2023/blob/main/scripts/processing/makeTadRanges.R
#              for how the source bedpe files were generated (Arrowhead,
#              SCALE-normalized, one file per timepoint).
#
#              Nested/hierarchical TAD calls (Arrowhead reports sub-TADs
#              inside larger domains) are flattened to only the outermost,
#              non-nested domain at each locus before computing sparsity,
#              so the null model isn't built on an arbitrarily small sub-TAD.
#
#              Sparsity is computed with clustAnalytics::internal_density(),
#              already a loopcity dependency (used for conductance elsewhere).
#              Given a graph and a vertex grouping, internal_density() returns
#              observed edges / all possible pairs within each group — i.e.
#              exactly the density of the naive-null complete graph's actual
#              fill, per TAD. The null's own density is 1 by construction, so
#              1 - internal_density() is directly "how much sparser."
#
# Input:       data/processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic
#              data/processed/hic/5kbLoops_0.txt (K562 loop calls)
#              data/processed/hic/tads/SCALE_0.bedpe (Arrowhead TAD calls,
#                0h condition — copy from MEGA2023 project data if not
#                already present in this project's data dir)
# Output:      figures/tad_null_model_sparsity.pdf
#              data/processed/tad_sparsity.csv
#              prints manuscript sentence to console
# ---------------------------------------------------------------------------


# Parameters ------------------------------------------------------------

project_dir <- "/work/users/j/p/jpflores/projects/loopcity"
data_dir    <- file.path(project_dir, "data")
output_dir  <- file.path(project_dir, "figures")

hic_file   <- file.path(data_dir, "processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic")
loops_file <- file.path(data_dir, "processed/hic/5kbLoops_0.txt")
tad_file   <- file.path(data_dir, "processed/hic/tads/SCALE_0.bedpe")

## Communities cache from the same pipeline used elsewhere — avoids
## re-running mergeAnchors/connectLoopAnchors/scoreInteractions/assignCommunities
communities_rds <- file.path(data_dir, "processed/loopcity_communities.rds")

hic_resolution <- 10e3
hic_norm       <- "NONE"
pixel_overlap  <- 1
overlap_dist   <- 1e6
leiden_res     <- 0.5


# Libraries ---------------------------------------------------------------

library(clustAnalytics)
library(data.table)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(glue)
library(gridExtra)
library(igraph)
library(InteractionSet)
library(loopcity)
library(mariner)


# Load / run loopcity pipeline ---------------------------------------------

if (file.exists(communities_rds)) {
  message("Loading cached communities from: ", communities_rds)
  communities <- readRDS(communities_rds)
} else {
  merged <- data.table::fread(loops_file) |>
    mariner::assignToBins(binSize = 10e3) |>
    mergeAnchors(pixelOverlap = pixel_overlap)
  GenomeInfoDb::seqlevelsStyle(merged) <- "ENSEMBL"
  connected   <- connectLoopAnchors(merged, overlapDist = overlap_dist)
  scored      <- scoreInteractions(connected,
                                   hicFile    = hic_file,
                                   resolution = hic_resolution,
                                   norm       = hic_norm)
  communities <- assignCommunities(InteractionSet::interactions(scored),
                                   leidenResolution = leiden_res)
  saveRDS(communities, communities_rds)
}


# Load TAD calls and flatten nested domains --------------------------------

## Arrowhead bedpe: V1:V3 = anchor1 (chr,start,end), V4:V6 = anchor2,
## V12 = score. A TAD's genomic span is anchor1$start to anchor2$end.
tad_raw <- data.table::fread(tad_file)
tad_gr <- GenomicRanges::GRanges(
  seqnames = tad_raw$V1,
  ranges   = IRanges::IRanges(start = tad_raw$V2 + 1, end = tad_raw$V6),
  score    = tad_raw$V12
)
GenomeInfoDb::seqlevelsStyle(tad_gr) <- "ENSEMBL"

## Keep only outermost domains: drop any TAD fully nested within a larger one.
## Iterate largest-to-smallest, keeping a TAD only if it isn't `within` an
## already-kept (necessarily larger) domain.
flattenNestedTads <- function(gr) {
  gr <- gr[order(-GenomicRanges::width(gr))]
  keep <- logical(length(gr))
  kept_gr <- GenomicRanges::GRanges()
  for (i in seq_along(gr)) {
    nested <- length(kept_gr) > 0 &&
      length(IRanges::subsetByOverlaps(gr[i], kept_gr, type = "within")) > 0
    if (!nested) {
      keep[i] <- TRUE
      kept_gr <- c(kept_gr, gr[i])
    }
  }
  gr[keep]
}

top_tads <- flattenNestedTads(tad_gr)
top_tads$tadId <- seq_along(top_tads)

message(glue(
  "{length(tad_gr)} raw TAD calls flattened to ",
  "{length(top_tads)} outermost (non-nested) domains"
))


# Build genome-wide anchor graph from final loopcity network ---------------

## Nodes = every merged anchor (regardless of degree); edges = interactions
## retained by assignCommunities (i.e. the network *after* community
## detection / pruning — the graph the reviewer wants compared to the null).
all_anchors <- InteractionSet::regions(communities)
anc <- InteractionSet::anchorIds(communities)

edge_df <- data.frame(from = anc$first, to = anc$second)
g <- igraph::graph_from_data_frame(
  edge_df,
  directed = FALSE,
  vertices = data.frame(name = seq_along(all_anchors))
)

## Assign each anchor to a flattened top-level TAD (NA if it falls outside
## all TAD calls — those anchors are excluded from the sparsity comparison)
overlaps <- GenomicRanges::findOverlaps(all_anchors, top_tads)
tad_membership <- rep(NA_integer_, length(all_anchors))
tad_membership[S4Vectors::queryHits(overlaps)] <-
  top_tads$tadId[S4Vectors::subjectHits(overlaps)]

## Restrict to anchors assigned to exactly one TAD
keep_vertices <- which(!is.na(tad_membership))
g_sub <- igraph::induced_subgraph(g, keep_vertices)
com_sub <- tad_membership[keep_vertices]


# Compute sparsity per TAD --------------------------------------------------

## internal_density = observed edges / all possible pairs within the group.
## The naive null (fully connected graph) has density 1 by construction, so
## 1 - density is exactly "how much sparser" loopcity's network is.
tad_density <- clustAnalytics::internal_density(g_sub, com_sub)
tad_edges   <- clustAnalytics::edges_inside(g_sub, com_sub)
tad_sizes   <- as.integer(table(factor(com_sub, levels = sort(unique(com_sub)))))

tad_sparsity <- data.frame(
  tadId          = sort(unique(com_sub)),
  n_anchors      = tad_sizes,
  possible_edges = choose(tad_sizes, 2),
  observed_edges = tad_edges,
  density        = tad_density
) |>
  dplyr::filter(n_anchors >= 2)

tad_sparsity$sparsity <- 1 - tad_sparsity$density

data.table::fwrite(tad_sparsity,
                   file.path(data_dir, "processed/tad_sparsity.csv"))


# Summary stats + manuscript sentence ---------------------------------------

med_density  <- median(tad_sparsity$density, na.rm = TRUE)
med_sparsity <- 1 - med_density

message(glue(
  "Median density across {nrow(tad_sparsity)} TADs: ",
  "{round(med_density, 4)} ",
  "(naive fully-connected null = 1)"
))

cat(glue(
  "Compared to a naive null model in which all merged loop anchors within ",
  "a TAD form a fully connected graph, loopcity's community-detected ",
  "network retained a median of only {round(med_density * 100, 2)}% of ",
  "possible anchor-anchor edges per TAD (n = {nrow(tad_sparsity)} TADs), ",
  "a {round(med_sparsity * 100, 1)}% reduction in edge density."
), "\n\n")


# Visualization -------------------------------------------------------------

## Panel 1: observed edges vs. naive null (choose(k,2)) as a function of
## TAD size — null grows quadratically, loopcity's kept edges grow far slower
p1 <- ggplot(tad_sparsity, aes(x = n_anchors)) +
  geom_function(fun = function(k) choose(k, 2),
                aes(color = "Naive null (fully connected)"),
                linewidth = 0.8) +
  geom_point(aes(y = observed_edges, color = "loopcity (post community detection)"),
             alpha = 0.5, size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c(
    "Naive null (fully connected)" = "grey60",
    "loopcity (post community detection)" = "deepskyblue3"
  )) +
  labs(x = "Merged anchors per TAD", y = "Edges",
       color = NULL,
       title = "Observed edges vs. naive fully-connected null, per TAD") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

## Panel 2: distribution of per-TAD density vs. the null's density of 1
p2 <- ggplot(tad_sparsity, aes(x = density)) +
  geom_histogram(bins = 40, fill = "deepskyblue3", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1, y = Inf, label = "naive null density = 1",
           hjust = 1.05, vjust = 1.5, size = 3, color = "grey40") +
  labs(x = "Internal density (observed edges / possible edges)",
       y = "Number of TADs",
       title = "Per-TAD network density after loopcity community detection") +
  theme_bw(base_size = 11)

pdf(file.path(output_dir, "tad_null_model_sparsity.pdf"), width = 7, height = 8)
gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()


# Session info --------------------------------------------------------------

sessionInfo()
