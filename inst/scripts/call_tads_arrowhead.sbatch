#!/bin/bash
#SBATCH -J arrowhead_tads_0h
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 240
#SBATCH --mem=32g
#SBATCH -o .batch_logs/arrowhead_tads_0h.%j.out
#SBATCH -e .batch_logs/arrowhead_tads_0h.%j.err
#-----------------------------------------------------------------------------
# Calls Arrowhead TADs (SCALE-normalized, 10kb) on the K562 MEGA 0h .hic file,
# matching the resolution used by inst/scripts/tad_null_model_sparsity.R.
# Output is renamed to SCALE_0.bedpe, the path that script expects.
#
# Run from the project root:
#   mkdir -p .batch_logs
#   sbatch inst/scripts/call_tads_arrowhead.sbatch
#-----------------------------------------------------------------------------
set -euo pipefail

JUICER_TOOLS=/proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_2.13.07.jar
HIC_FILE=/work/users/j/p/jpflores/projects/loopcity/data/processed/hic/GSE214123_MEGA_K562_WT_0_inter.hic
OUT_DIR=/work/users/j/p/jpflores/projects/loopcity/data/processed/hic/tads

mkdir -p "$OUT_DIR"

java -jar "$JUICER_TOOLS" arrowhead -r 10000 -k SCALE "$HIC_FILE" "$OUT_DIR"

mv "$OUT_DIR/10000_blocks.bedpe" "$OUT_DIR/SCALE_0.bedpe"

echo "Done. TAD calls written to $OUT_DIR/SCALE_0.bedpe"
