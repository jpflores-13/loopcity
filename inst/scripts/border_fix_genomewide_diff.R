# Border-adjustment bug fix: genome-wide before/after diff ----------------
# Author:      JP Flores
# Date:        2026-07-13
# Project:     loopcity
# Description: Quantifies the real-data impact of the assignCommunities()
#              border-adjustment fix (commits fdc1731, 527eb78, ce26123),
#              addressing PI question "are we sure this doesn't change the
#              results?" The chr22 toy-data test (188 loops) showed the fix
#              is not a no-op (1 previously-orphaned loop correctly rescued),
#              but that dataset is too small to know the real magnitude.
#
#              Isolates ONLY the border-adjustment logic: mergeAnchors,
#              connectLoopAnchors, and scoreInteractions are run once via the
#              current package source (identical for both arms), and the
#              RNG seed is pinned immediately before each assignCommunities()
#              call so Leiden clustering draws the same random sequence.
#              Only .compareNodeToComm / assignCommunities differ between
#              the "old" (pre-fix, commit b252df2) and "new" (current) runs.
#
#              Also flags whether any changed loop falls inside the
#              chr1:110-115Mb region used for Figure 1b/c in
#              survey_communities.R, since that's the one figure already
#              built from this pipeline's output.
#
# Input:       data/processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic
#              data/processed/hic/5kbLoops_0.txt (K562 loop calls)
# Output:      data/processed/border_fix_diff.csv (all loops whose
#                loopCommunity assignment changed)
#              prints summary + manuscript sentence to console
# ---------------------------------------------------------------------------


# Parameters ----------------------------------------------------------------

## Path to the loopcity package SOURCE checkout (a git clone of
## jpflores-13/loopcity with full history) — used only to pull the pre-fix
## versions of assignCommunities.R/utils.R via `git show`. This is likely a
## different path than `project_dir` below if data and package source live
## in separate directories on the cluster.
loopcity_repo_dir <- "/work/users/j/p/jpflores/projects/loopcity"

project_dir <- "/work/users/j/p/jpflores/projects/loopcity"
data_dir    <- file.path(project_dir, "data")

hic_file   <- file.path(data_dir, "processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic")
loops_file <- file.path(data_dir, "processed/hic/5kbLoops_0.txt")

hic_resolution <- 10e3
hic_norm       <- "NONE"
pixel_overlap  <- 1
overlap_dist   <- 1e6
leiden_res     <- 0.5   ## matches survey_communities.R / benchmark_timing.R

pre_fix_commit <- "b252df2"   ## last commit before the border-adjustment fixes

## Figure 1b/c region (survey_communities.R) — flagged separately below
survey_chrom      <- "1"
survey_chromstart <- 110e6
survey_chromend   <- 115e6


# Libraries -------------------------------------------------------------

suppressPackageStartupMessages({
  devtools::load_all(loopcity_repo_dir, quiet = TRUE)
  library(data.table)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(glue)
  library(InteractionSet)
  library(mariner)
})


# Run shared pipeline steps ONCE (identical for both arms) -----------------

message("Running mergeAnchors -> connectLoopAnchors -> scoreInteractions ...")

merged <- data.table::fread(loops_file) |>
  mariner::assignToBins(binSize = 10e3) |>
  mergeAnchors(pixelOverlap = pixel_overlap)
GenomeInfoDb::seqlevelsStyle(merged) <- "ENSEMBL"

connected <- connectLoopAnchors(merged, overlapDist = overlap_dist)
scored    <- scoreInteractions(connected,
                               hicFile    = hic_file,
                               resolution = hic_resolution,
                               norm       = hic_norm) |>
  InteractionSet::interactions()


# Run 1: current (fixed) assignCommunities -----------------------------------

fixed_assignCommunities <- assignCommunities

set.seed(1)
new_result <- fixed_assignCommunities(scored, leidenResolution = leiden_res)
message(glue("NEW (fixed) run: {length(new_result)} loops"))


# Run 2: pre-fix assignCommunities, pulled from git and sourced to shadow ----

tmp_dir <- tempdir()
old_utils_path <- file.path(tmp_dir, "old_utils.R")
old_assign_path <- file.path(tmp_dir, "old_assignCommunities.R")

get_old_file <- function(relpath, out_path) {
  ok <- system2("git",
                args = c("-C", shQuote(loopcity_repo_dir), "show",
                          glue("{pre_fix_commit}:{relpath}")),
                stdout = out_path)
  if (!isTRUE(ok == 0) || file.size(out_path) == 0) {
    rlang::abort(glue(
      "Could not retrieve {relpath} at {pre_fix_commit} from ",
      "{loopcity_repo_dir}. Make sure that directory is a git checkout ",
      "with this commit (try `git fetch origin` there)."
    ))
  }
}
get_old_file("R/utils.R", old_utils_path)
get_old_file("R/assignCommunities.R", old_assign_path)

source(old_utils_path)
source(old_assign_path)
stopifnot(!identical(environment(assignCommunities),
                     environment(fixed_assignCommunities)))

set.seed(1)
old_result <- assignCommunities(scored, leidenResolution = leiden_res)
message(glue("OLD (pre-fix) run: {length(old_result)} loops"))


# Compare ---------------------------------------------------------------

mk_key <- function(x) paste(
  GenomeInfoDb::seqnames(InteractionSet::anchors(x, "first")),
  start(InteractionSet::anchors(x, "first")),
  start(InteractionSet::anchors(x, "second"))
)

old_df <- data.frame(
  key = mk_key(old_result),
  chr = as.character(GenomeInfoDb::seqnames(InteractionSet::anchors(old_result, "first"))),
  start1 = start(InteractionSet::anchors(old_result, "first")),
  start2 = start(InteractionSet::anchors(old_result, "second")),
  loopCommunity_old = I(lapply(old_result$loopCommunity, sort))
)
new_df <- data.frame(
  key = mk_key(new_result),
  loopCommunity_new = I(lapply(new_result$loopCommunity, sort))
)

merged_df <- merge(old_df, new_df, by = "key")
changed <- mapply(function(a, b) !identical(a, b),
                  merged_df$loopCommunity_old, merged_df$loopCommunity_new)

diff_df <- merged_df[changed, c("chr", "start1", "start2",
                                "loopCommunity_old", "loopCommunity_new")]
diff_df$loopCommunity_old <- vapply(diff_df$loopCommunity_old,
                                    paste, character(1), collapse = ",")
diff_df$loopCommunity_new <- vapply(diff_df$loopCommunity_new,
                                    paste, character(1), collapse = ",")
diff_df$rescued <- diff_df$loopCommunity_old == ""

## Flag whether any changed loop falls in the Figure 1b/c survey region
in_survey_region <- diff_df$chr == survey_chrom &
  diff_df$start1 >= survey_chromstart & diff_df$start2 <= survey_chromend
diff_df$in_fig1bc_region <- in_survey_region

data.table::fwrite(diff_df, file.path(data_dir, "processed/border_fix_diff.csv"))


# Summary + manuscript sentence -------------------------------------------

n_total   <- nrow(merged_df)
n_changed <- sum(changed)
n_rescued <- sum(diff_df$rescued)
n_lost    <- sum(!diff_df$rescued & diff_df$loopCommunity_new == "")
n_in_fig  <- sum(diff_df$in_fig1bc_region)

cat(glue("\n=== Genome-wide border-adjustment fix impact ===\n"))
cat(glue("Total loops compared: {n_total}\n"))
cat(glue("Loops with changed loopCommunity assignment: {n_changed} ",
         "({round(100 * n_changed / n_total, 3)}%)\n"))
cat(glue("  - Rescued (previously unassigned, now assigned): {n_rescued}\n"))
cat(glue("  - Lost (previously assigned, now unassigned): {n_lost}\n"))
cat(glue("  - Changed loops overlapping Fig 1b/c region ",
         "(chr1:110-115Mb): {n_in_fig}\n\n"))

cat(glue(
  "During revision we identified and fixed a bug in assignCommunities()'s ",
  "border-adjustment step ... [see commits fdc1731/527eb78/ce26123]. ",
  "Comparing pre-fix and post-fix outputs on the full genome-wide K562 ",
  "loop set (identical inputs, seed pinned), {n_changed} of {n_total} ",
  "loops ({round(100 * n_changed / n_total, 3)}%) changed community ",
  "assignment, {n_rescued} of which were previously-orphaned border loops ",
  "now correctly assigned. {ifelse(n_in_fig == 0, ",
  "'None of these loops fall within the region shown in Figure 1b/c.', ",
  "glue('{n_in_fig} of these loops fall within the region shown in ",
  "Figure 1b/c and that panel should be regenerated.'))}"
), "\n\n")

sessionInfo()
