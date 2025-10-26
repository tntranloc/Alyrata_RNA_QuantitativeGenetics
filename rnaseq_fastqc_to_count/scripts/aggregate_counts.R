#!/usr/bin/env Rscript

# =========================
# Aggregate STAR GeneCounts
# =========================
# Usage:
#   Rscript scripts/aggregate_counts.R results/mapping/ counts_matrix.csv
#
# Each STAR output file must be named:
#   SAMPLE.ReadsPerGene.out.tab
# Columns:
#   1: Gene ID
#   2: unstranded counts
#   3: read strand 1 (htseq -s yes)
#   4: read strand 2 (htseq -s reverse)
#
# Pick the correct column depending on your strandedness (set below)
# =========================

suppressPackageStartupMessages({
  library(tidyverse)
})

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript aggregate_counts.R <input_dir> <output_file>")
}

indir  = args[1]
outfile = args[2]

# -------------------------
# User parameter
# -------------------------
# choose which column of STAR output to use:
# 2 = unstranded, 3 = first strand, 4 = second strand
# Adjust according to your library prep (check FastQC or metadata)
count_col = 4  # <- change this if needed

# -------------------------
# Load all count tables
# -------------------------
files = list.files(indir, pattern = "ReadsPerGene.out.tab$", full.names = TRUE)
if (length(files) == 0) stop("No STAR ReadsPerGene.out.tab files found in ", indir)

message("Found ", length(files), " count tables.")

count_list = lapply(files, function(f) {
  sample_id = strsplit(basename(f), "\\.")[[1]][1]
  dat = read.delim(f, header = FALSE)
  dat = dat[, c(1, count_col)]
  colnames(dat) = c("gene", sample_id)
  return(dat)
})

# -------------------------
# Merge into one matrix
# -------------------------
counts = reduce(count_list, full_join, by = "gene")
counts[is.na(counts)] = 0

# Remove STAR summary rows (those starting with "__")
counts = counts[!grepl("^__", counts$gene), ]

# -------------------------
# Save results
# -------------------------
write.csv(counts, outfile, row.names = FALSE, quote = FALSE)
message("✅ Counts written to: ", outfile)

# Optional: print summary
message("Matrix dimensions: ", nrow(counts), " genes × ", ncol(counts) - 1, " samples")
