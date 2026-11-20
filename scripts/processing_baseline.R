# ============================================================================
# Baseline preparation (no signal added): writes the inputs needed to run
# every CCC method on the original processed data, per dataset x strategy.
#
# The baseline is compared with each semi-simulated run in the metric step:
# non-simulated interactions that become significant only after signal
# injection, and share no gene with the simulated pair, are false positives.
#
# Outputs mirror the files written by processing_semisimulation_nbsampling.R,
# so the method scripts can be reused unchanged:
#   - baseline_counts: original processed counts (normalized downstream by
#     processing_normalization.R, like the inflated counts)
#   - baseline_cellmetadata: processed cell metadata plus a Celltype_updated
#     column. Methods infer interactions from Celltype (Sender/Receiver), so
#     Celltype_updated is only used for their diagnostic columns; here ALL
#     Sender cells are labelled "Sender_signalAdded" and ALL Receiver cells
#     "Receiver_signalAdded" so those diagnostics run without errors. The
#     diagnostic columns of baseline outputs are not used downstream.
#   - metadata_cpdbv5: CellPhoneDB-formatted metadata (Celltype only)
# ============================================================================

suppressMessages({
  library(dplyr)
  library(magrittr)
  library(jsonlite)
  source("scripts/helper_functions.R")
})

# Fail early with a clear error if a required input is missing
if (is.null(snakemake@input[["processed_counts"]]) | is.null(snakemake@input[["cellmetadata"]])) {
  stop("Argument_name needs to be specified, but is missing.n", call. = FALSE)
}

# ==============================================================================
# STEP 0: Snakemake I/O
# ==============================================================================

processed_counts_path = snakemake@input[["processed_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata"]]

baseline_counts_path = snakemake@output[["baseline_counts"]]
baseline_cellmetadata_path = snakemake@output[["baseline_cellmetadata"]]
metadata_cpdbv5_path = snakemake@output[["metadata_cpdbv5"]]

# ==============================================================================
# STEP 1: load data (same way as processing_semisimulation_nbsampling.R)
# ==============================================================================

counts = read.table(processed_counts_path)
cellmetadata = read_json(path = cellmetadata_path)
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

# sanity check: cell names of counts and metadata must correspond and be in the same order
stopifnot(colnames(counts) == cellmetadata$metadata$Cell_ID)

# ==============================================================================
# STEP 2: label all Sender/Receiver cells as "signalAdded" (diagnostics only)
# ==============================================================================

cellmetadata$metadata$Celltype_updated = cellmetadata$metadata$Celltype
cellmetadata$metadata$Celltype_updated[cellmetadata$metadata$Celltype == "Sender"] = "Sender_signalAdded"
cellmetadata$metadata$Celltype_updated[cellmetadata$metadata$Celltype == "Receiver"] = "Receiver_signalAdded"

# ==============================================================================
# STEP 3: save outputs (same formats as the semi-simulation step)
# ==============================================================================

write.table(counts, baseline_counts_path, sep = "\t")
write_json(cellmetadata, baseline_cellmetadata_path)

# CellPhoneDB strictly requires a specific metadata file format
mt = cellmetadata$metadata %>% select(c("Cell_ID", "Celltype")) %>% setNames(c("barcode_sample", "cell_type"))
rownames(mt) = mt$barcode_sample
write.table(mt, metadata_cpdbv5_path, sep = "\t")
