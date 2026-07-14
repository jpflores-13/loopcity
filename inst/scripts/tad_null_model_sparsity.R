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
## skip = 2: juicer_tools 2.x prefixes output with a "#chr1 x1 x2..." column
## header line and a "# juicer_tools version..." comment line, both of which
## confuse fread's header autodetection if left in.
tad_raw <- data.table::fread(tad_file, skip = 2, header = FALSE)
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

## clustAnalytics's community functions require contiguous 1:k labels (they
## index communities positionally); tadId has gaps because many TADs have no
## overlapping anchors. Remap to contiguous labels for the call, keeping
## tad_ids_present to translate back to genomic TAD identity in the output.
tad_ids_present <- sort(unique(com_sub))
com_contig <- match(com_sub, tad_ids_present)

## clustAnalytics::internal_density()/edges_inside() build a weighted copy of
## the graph internally but then (a package bug) compute the edge list from
## the original, unweighted graph if it has no "weight" attribute -- silently
## returning garbage values instead of erroring. Set weight explicitly so
## they never hit that path.
g_sub <- igraph::set_edge_attr(g_sub, "weight", value = 1)


# Compute sparsity per TAD --------------------------------------------------

## internal_density = observed edges / all possible pairs within the group.
## The naive null (fully connected graph) has density 1 by construction, so
## 1 - density is exactly "how much sparser" loopcity's network is.
tad_density <- clustAnalytics::internal_density(g_sub, com_contig)
tad_edges   <- clustAnalytics::edges_inside(g_sub, com_contig)
tad_sizes   <- as.integer(table(factor(com_contig, levels = seq_along(tad_ids_present))))

tad_sparsity <- data.frame(
  tadId          = tad_ids_present,
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

mean_density  <- mean(tad_sparsity$density, na.rm = TRUE)
mean_sparsity <- 1 - mean_density

message(glue(
  "Mean density across {nrow(tad_sparsity)} TADs: ",
  "{round(mean_density, 4)} ",
  "(naive fully-connected null = 1)"
))

cat(glue(
  "Compared to a naive null model in which all merged loop anchors within ",
  "a TAD form a fully connected graph (mean density = 1 by construction), ",
  "loopcity's community-detected network had a mean density of only ",
  "{round(mean_density, 2)} across {nrow(tad_sparsity)} TADs, a ",
  "{round(mean_sparsity * 100, 1)}% reduction in edge density relative to ",
  "the naive null."
), "\n\n")


# Visualization -------------------------------------------------------------

loopcity_color <- "#2a78d6"  # blue -- loopcity (observed)

## Density by TAD-size bin vs. the null's density of 1. A raw
## histogram of density is misleading here: small TADs have few possible
## anchor pairs (choose(k,2)), so density can only take a handful of discrete
## values (e.g. exactly 0 or 1 for a 2-anchor TAD) -- producing spiky,
## multimodal-looking bars that are a combinatorial artifact, not biology.
## Binning by TAD size and boxplotting shows sparsity holds across size
## classes, not just among small, discretization-prone TADs.
tad_sparsity$size_bin <- cut(
  tad_sparsity$n_anchors,
  breaks = c(1, 3, 6, 10, Inf),
  labels = c("2-3", "4-6", "7-10", "11+")
)

## Direct-label each box with its median so the headline numbers are legible
## without having to visually read quartiles off the plot.
box_medians <- tad_sparsity |>
  dplyr::group_by(size_bin) |>
  dplyr::summarize(med = median(density), .groups = "drop")

p2 <- ggplot(tad_sparsity, aes(x = size_bin, y = density)) +
  geom_boxplot(fill = loopcity_color, alpha = 0.6, outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  annotate("text", x = Inf, y = 1, label = "naive null density = 1",
           hjust = 1.05, vjust = 1.5, size = 3, color = "grey40") +
  geom_text(data = box_medians, aes(x = size_bin, y = med, label = sprintf("%.2f", med)),
            vjust = -1, size = 3.2, fontface = "bold", color = "grey20",
            inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(x = "Merged anchors per TAD", y = "Internal density (observed / possible edges)",
       title = "Network density by TAD size") +
  theme_bw(base_size = 11)

pdf(file.path(output_dir, "tad_null_model_sparsity.pdf"), width = 6, height = 4.5)
print(p2)
dev.off()


# Session info --------------------------------------------------------------

sessionInfo()
