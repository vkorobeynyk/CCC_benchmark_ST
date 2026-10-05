# ============================================================================
# Semi-simulation step: picks one ligand-receptor pair (indexLR_toSample),
# inflates its ligand expression in a subset of Sender cells and its receptor
# expression in a subset of Receiver cells (via semi_simulate(), negative-
# binomial sampling), and writes out the inflated counts plus updated cell
# metadata for downstream CCC methods.
#
# The number of Sender/Receiver cells that receive inflated expression is
# controlled by PCE_Sender / PCE_Receiver, which multiply the dataset's
# average percentage of cells expressing L/R genes (see config.yaml for the
# naming caveat).
# ============================================================================

suppressMessages({
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(ggpubr)
  library(sf)
  library(jsonlite)
  source("scripts/helper_functions.R")
})

# Fail early with a clear error if a required input is missing
if (is.null(snakemake@input[["processed_counts"]]) | is.null(snakemake@input[["genemetadata"]]) | is.null(snakemake@params[["LR_database"]]) |
    is.null(snakemake@input[["cellmetadata"]])) {
  stop("Argument_name needs to be specified, but is missing.n", call. = FALSE)
}

# ==============================================================================
# STEP 0: Snakemake I/O and params
# ==============================================================================

# input files
processed_counts_path = snakemake@input[["processed_counts"]]
genemetadata_path = snakemake@input[["genemetadata"]]
cellmetadata_path = snakemake@input[["cellmetadata"]]

# output files
inflated_counts_path = snakemake@output[["inflated_counts"]]
simulated_cellmetadata_path = snakemake@output[["simulated_cellmetadata"]]
simulated_interactions_path = snakemake@output[["simulated_interactions"]]
metadata_cpdbv5_path = snakemake@output[["metadata_cpdbv5"]]
plot_neighbors_path = snakemake@output[["plot_neighbors"]]

# params
indexLR_toSample = snakemake@wildcards[["indexLR_toSample"]] %>% as.integer
PCE_Sender = snakemake@wildcards[["PCE_Sender"]] %>% as.double
PCE_Receiver = snakemake@wildcards[["PCE_Receiver"]] %>% as.double
LR_database_path = snakemake@params[["LR_database"]]

# ==============================================================================
# STEP 1: load data
# ==============================================================================

counts = read.table(processed_counts_path)
genemetadata = readRDS(genemetadata_path)
means_perCT = genemetadata$mean
cellmetadata = read_json(path = cellmetadata_path)

# ---------------------------------------------------------------------------
# Manual interactive testing snippet (not executed): quickly load a single
# dataset's processed files in an R session, e.g. for debugging outside Snakemake.
# ---------------------------------------------------------------------------
'
counts = read.table("data/processed/spatialScattering/CosMx_HFC/processed_counts_CosMx_HFC.tsv")
genemetadata = readRDS("data/processed/spatialScattering/CosMx_HFC/genemetadata_CosMx_HFC.RDS")
cellmetadata = read_json("data/processed/spatialScattering/CosMx_HFC/cellmetadata_CosMx_HFC.json")
means_perCT = genemetadata$mean
LRdb = read.table("data/LR_database.tsv", header = T)
'

# transform the JSON list into individual dataframes
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

# NOTE: cellmetadata contains all celltypes and all cells. At the end of this
# script, cellmetadata is saved as the metadata that the CCC methods will take as input.

# sanity check: cell names of counts and metadata must correspond and be in the same order
stopifnot(colnames(counts) == cellmetadata$metadata$Cell_ID)

message(paste("Shape of counts object:", str_flatten(dim(counts), " ")))

# ==============================================================================
# STEP 2: prepare the L/R database and gene mean table for this run
# ==============================================================================

# drop genes with mean == 0 in either celltype to be simulated, since they
# can't be used to parameterize the negative-binomial sampling below
gene_index = which(means_perCT[, "Sender"] == 0 | means_perCT[, "Receiver"] == 0)
if (length(gene_index) > 0) {
  means_perCT = means_perCT[-gene_index, ]
}

LRdb = read.table(LR_database_path, header = T)

# keep only LR pairs whose genes are present in the count data and have a non-zero mean per celltype
n = LRdb %>%
  apply(., 1, function(row) {
    str_split(row, "_") %>%
      lapply(., function(x) {(x %in% rownames(counts)) & (x %in% means_perCT$gene_names)}) %>%
      unlist %>%
      all
  })
LRdb = LRdb[n, ]

message(paste("After filtering LR database,", nrow(LRdb), "LR pairs show expression in at least 10 cells"))

# ==============================================================================
# STEP 3: select the ligand-receptor pair to inflate expression for
#
# Selection is driven by indexLR_toSample. Half of the sampled LR pairs have
# a single-subunit receptor, half have a 2-subunit receptor.
# ==============================================================================

simulated_interactions_lst = list()
comb_CT = "Sender_Receiver" # pre-sampled LR pair (and its subunit genes) is stored under this key

config = yaml::read_yaml("config.yaml")

# split the LR database by number of receptor subunits
LR_subset_1_under = LRdb[str_count(LRdb$ligand_receptor, "_") == 1, ]
LR_subset_2_under = LRdb[str_count(LRdb$ligand_receptor, "_") == 2, ]

idx = match(indexLR_toSample, config$indexLR_toSample)
#halfway = length(config$indexLR_toSample) / 2
halfway = length(config$indexLR_toSample) / 2


if (idx <= halfway) {
  # FIRST HALF: single-subunit L/R pairs
  LR_sample = LR_subset_1_under[indexLR_toSample, ]
  
} else {
  # SECOND HALF: 2-subunit L/R pairs
  local_idx = idx - halfway
  
  # saveguard against indexLR_toSample being larger than the available 2-subunit pairs
  stopifnot((indexLR_toSample / 2) <= nrow(LR_subset_2_under))
  
  LR_sample = LR_subset_2_under[local_idx, ]
}

simulated_interactions_lst[[comb_CT]]$ligand = LR_sample$ligand
simulated_interactions_lst[[comb_CT]]$receptor = LR_sample$receptor

# ==============================================================================
# STEP 4: inflate ligand/receptor expression (semi_simulate)
# ==============================================================================

# Scale the dataset's average percentage of cells expressing L/R genes by
# PCE_Sender / PCE_Receiver to get the target fraction of Sender/Receiver
# cells that should express the ligand/receptor after inflation.
fraction_cells_expressingR = PCE_Receiver * cellmetadata$average_percentageCells_expressingLR
fraction_cells_expressingL = PCE_Sender * cellmetadata$average_percentageCells_expressingLR

# if these fractions exceed 1, we'd be trying to sample more cells than are available -
# decrease PCE_Receiver or PCE_Sender in config.yaml
stopifnot(fraction_cells_expressingR < 1 & fraction_cells_expressingL < 1)

# subunits are simulated too, since the LR database nomenclature can be L_R1_R2 etc.
semi_simulation_out = semi_simulate(counts = counts, simulated_interactions_lst = simulated_interactions_lst, genemetadata = genemetadata,
                                    metadata = cellmetadata$metadata,
                                    fraction_cells_expressingR = fraction_cells_expressingR,
                                    fraction_cells_expressingL = fraction_cells_expressingL,
                                    df_neighbors = cellmetadata$neighbor_cells)

# ==============================================================================
# STEP 5: label signal-added cells and produce a diagnostic plot
# ==============================================================================

cellmetadata$metadata$Celltype_updated = cellmetadata$metadata$Celltype
cellmetadata$metadata$Celltype_updated[cellmetadata$metadata$Cell_ID %in% semi_simulation_out$cell_info$Sender_Receiver$cells_signalAdded$Sender] = "Sender_signalAdded"
cellmetadata$metadata$Celltype_updated[cellmetadata$metadata$Cell_ID %in% semi_simulation_out$cell_info$Sender_Receiver$cells_signalAdded$Receiver] = "Receiver_signalAdded"

p = ggplot(cellmetadata$metadata, aes(x = x, y = y, color = Celltype_updated, size = Celltype_updated)) +
  geom_point() +
  xlab("x") +
  ylab("y") +
  scale_color_manual(values = c("#FFCCFF", "#990099", "#CCCCFF", "#0000FF", "#FFCC99")) +
  scale_size_manual(values = c(1.5, 3, 1.5, 3, 0.5)) +
  theme_bw()

ggsave(filename = plot_neighbors_path, plot = p, width = 200, height = 150, units = "mm")

# ==============================================================================
# STEP 6: save outputs
# ==============================================================================

write.table(semi_simulation_out$counts_inflated, inflated_counts_path, sep = "\t")
saveRDS(simulated_interactions_lst, simulated_interactions_path)

write_json(cellmetadata, simulated_cellmetadata_path)

# CellPhoneDB strictly requires a specific metadata file format - create it here
mt = cellmetadata$metadata %>% select(c("Cell_ID", "Celltype")) %>% setNames(c("barcode_sample", "cell_type"))
rownames(mt) = mt$barcode_sample
write.table(mt, metadata_cpdbv5_path, sep = "\t")