# ============================================================================
# Normalizes the semi-simulated ("inflated") counts using scater's
# logNormCounts, so every CCC method downstream consumes the same normalized
# expression matrix.
# ============================================================================

suppressMessages({
  library(scater)
  library(SpatialExperiment)
  library(dplyr)
})

# Fail early with a clear error if a required input/output is missing
if (is.null(snakemake@input[["inflated_counts"]]) | is.null(snakemake@output[["normalized_counts"]])) {
  stop("Argument_name needs to be specified, but is missing.n", call. = FALSE)
}

# ==============================================================================
# STEP 0: Snakemake I/O
# ==============================================================================

# input files
inflated_counts_path = snakemake@input[["inflated_counts"]]

# output files
normalized_counts_path = snakemake@output[["normalized_counts"]] 

# ==============================================================================
# STEP 1: load the inflated counts
# ==============================================================================

inflated_counts = read.table(inflated_counts_path)

# ==============================================================================
# STEP 2: log-normalize (scater::logNormCounts)
# ==============================================================================

normalized_inflated_counts = SingleCellExperiment(list(counts = inflated_counts)) %>%
  logNormCounts %>%
  assay(., "logcounts")

# ==============================================================================
# STEP 3: save output
# ==============================================================================

write.table(normalized_inflated_counts, normalized_counts_path, sep = "\t")