# ============================================================================
# First step of the spatial CCC benchmark pipeline, "colocalization" setting:
# Sender and Receiver cells are drawn from a manually curated set of cells
# (selected manually via the shiny app in
# /CCC_benchmark_ST/shiny_forSelectingCells) that sit close together in a
# specific spatial region of interest, then randomly split in half into
# Sender and Receiver. 
#
# For a given spatial dataset, this script:
#   1. Runs generic QC (general_ST_QCpipeline_R) without filtering, only to
#      produce diagnostic plots.
#   2. Loads the manually selected cells for this dataset and splits them
#      into Sender / Receiver.
#   3. Estimates the dataset's average percentage of cells expressing L/R genes
#      (used later to calibrate the semi-simulation).
#   4. Pairs each Sender cell with its nearest available Receiver cell
#      (find_neighboring_spots); since cells were already chosen to colocalize,
#      the distance filtering bounds are left wide open ([0, Inf]).
#   5. Estimates per-gene mean/dispersion (edgeR) to be used for the
#      negative-binomial semi-simulation.
#   6. Writes out processed counts, gene metadata, cell metadata (as JSON),
#      and several diagnostic plots.
# ============================================================================

suppressMessages({
  library(dplyr)
  library(magrittr)
  library(edgeR)
  library(sf)
  library(stringr)
  library(ggspavis)
  library(scater)
  library(SingleCellExperiment)
  library(jsonlite)
  library(ggExtra)
  library(ggforce)
  source("scripts/helper_functions.R")
  source("scripts/general_spatial_QCpipeline_R_function.R")
})

# Fail early with a clear error if a required input is missing
if (is.null(snakemake@input[["counts"]]) | is.null(snakemake@input[["metadata"]])) {
  stop("Argument_name needs to be specified, but is missing.n", call. = FALSE)
}

# ==============================================================================
# STEP 0: Snakemake I/O and params
# ==============================================================================

# input files
counts_path = snakemake@input[["counts"]]
metadata_path = snakemake@input[["metadata"]]

# output files
processed_counts_path = snakemake@output[["processed_counts"]]
genemetadata_path = snakemake@output[["genemetadata"]]
cellmetadata_path = snakemake@output[["cellmetadata"]]
plot_allneighbors_path = snakemake@output[["plot_allneighbors"]]
plot_radius = snakemake@output[["plot_radius"]]
plot_SenderReceiverDistance = snakemake@output[["plot_SenderReceiverDistance"]]
diagnostic_plots_path = snakemake@output[["diagnostic_plots"]]

# params
LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)
dataset = snakemake@params["dataset"] %>% unlist

# ==============================================================================
# STEP 1: load data and run QC diagnostics only (no filtering)
# ==============================================================================

counts = read.table(counts_path, row.names = 1)
metadata = read.table(metadata_path, row.names = 1)

# run the QC pipeline without any filtering - we only want the diagnostic plots here
x = general_ST_QCpipeline_R(counts, to_filter = FALSE)

pdf(diagnostic_plots_path)
x[-1]
dev.off()

# ---------------------------------------------------------------------------
# Manual interactive testing snippet (not executed): quickly load a single
# dataset's raw files in an R session, e.g. for debugging outside Snakemake.
# ---------------------------------------------------------------------------
'
counts = read.table("data/MERFISH_mColon/counts_MERFISH_mColon.tsv", row.names = 1)
metadata = read.table("data/MERFISH_mColon/metadata_MERFISH_mColon.tsv",row.names = 1)
LR_database = read.table("data/LR_database.tsv", row.names = 1)
counts = read.table("data/Visium_HD_HPC/counts_Visium_HD_HPC.tsv", row.names = 1)
metadata = read.table("data/Visium_HD_HPC/metadata_Visium_HD_HPC.tsv",row.names = 1)
LR_database = read.table("data/LR_database.tsv", row.names = 1)
counts = read.table("data/CosMx_HFC/counts_CosMx_HFC.tsv", row.names = 1)
metadata = read.table("data/CosMx_HFC/metadata_CosMx_HFC.tsv",row.names = 1)
LR_database = read.table("data/LR_database.tsv", row.names = 1)
'

# remove genes with 0 counts and keep genes expressed in at least 10 cells
counts = counts[rowSums(counts) != 0 & rowSums(counts != 0) > 10, ]

# ==============================================================================
# STEP 2: load the manually selected, colocalized cells and split them into
# Sender / Receiver
# ==============================================================================

set.seed(1)
if (dataset == "Visium_HD_HPC") {cells = read.csv("data/Visium_HD_HPC/DG_cells.csv") %>% pull(Cell_ID)
} else if (dataset == "MERFISH_mColon") {cells = read.csv("data/MERFISH_mColon/selected_cells.csv") %>% pull(Cell_ID)
} else if (dataset == "CosMx_HFC") {cells = read.csv("data/CosMx_HFC/selected_cells.csv") %>% pull(Cell_ID)
}

# split the manually selected cells in half: Sender and Receiver
sender_cells = sample(cells, length(cells) * 0.5)
receiver_cells = setdiff(cells, sender_cells)

metadata[, "Celltype"] = "Other"
metadata$Celltype[metadata$Cell_ID %in% sender_cells] = "Sender"
metadata$Celltype[metadata$Cell_ID %in% receiver_cells] = "Receiver"

# ==============================================================================
# STEP 3: estimate the dataset's average % of cells expressing L/R genes
#
# This value is the starting point for calibrating the simulated fraction of
# cells/spots expressing L/R genes later on
# ==============================================================================

LRgenes = c(LR_database$ligand, LR_database$receptor) %>%
  str_split("_") %>%
  unlist %>%
  unique()

average_percentageCells_expressingLR = apply(counts[which(rownames(counts) %in% LRgenes), ], 1, function(x) {
  sum(x > 0) / length(x)
}) %>% mean

# ==============================================================================
# STEP 4: pair each Sender cell with its nearest available Receiver cell
#
# Since Sender/Receiver cells were already hand-picked to colocalize, the
# post-filtering distance bounds are left wide open ([0, Inf]) - we don't
# need to additionally filter out far-away pairs here the way the scattered-
# signal setting does.
# ==============================================================================

if (dataset == "Visium_HD_HPC") {
  distance = 250
  distance_post_filtering_lower = 0
  distance_post_filtering_upper = Inf
} else if (dataset == "MERFISH_mColon") {
  distance = 50
  distance_post_filtering_lower = 0
  distance_post_filtering_upper = Inf
} else if (dataset == "CosMx_HFC") {
  distance = 250
  distance_post_filtering_lower = 0
  distance_post_filtering_upper = Inf
}

neighbors_info = find_neighboring_spots(spatial_coords = metadata %>% select(c("x", "y")),
                                        ligand_spots = metadata %>% filter(Celltype == "Sender") %>% select(Cell_ID) %>% unlist %>% unname,
                                        receptor_spots = metadata %>% filter(Celltype == "Receiver") %>% select(Cell_ID) %>% unlist %>% unname,
                                        distance_threshold_dataset = distance,
                                        distance_post_filtering_lower = distance_post_filtering_lower,
                                        distance_post_filtering_upper = distance_post_filtering_upper)

table(neighbors_info$metadata$Celltype)

# subset counts to only the cells kept in metadata
counts = counts[, metadata$Cell_ID]

# diagnostic plot: distance from each Sender to its closest Receiver
df = data.frame(d = neighbors_info$distance_Sender_toClosest_Receiver)

p = ggplot(df, aes(x = seq_along(d), y = d)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(x = "Index of Sender cell", y = "Distance", title = "Distance to Closest Receiver")

ggsave(filename = plot_SenderReceiverDistance, plot = p, width = 200, height = 150, units = "mm")

# ==============================================================================
# STEP 5: finalize Sender/Receiver/Other labels and visualize the layout
# ==============================================================================

# update metadata according to the pairing found by find_neighboring_spots()
metadata = neighbors_info$metadata
neighbor_cells = neighbors_info$neighbors

# In order to have the same number of cells as other semi-simulation settings,
# optionally subsample to a fixed N_cells
#N_cells = 100
#neighbor_cells = neighbor_cells[,sample(neighbor_cells, N_cells) %>% names]

# diagnostic plot: simple Sender/Receiver/Other spatial scatter
p = ggplot() +
  # draw "Other" cells first, small and faded, so they read as background
  geom_point(data = metadata %>% filter(Celltype == "Other"),
             aes(x = x, y = y), color = "grey75", size = 0.6, alpha = 0.6) +
  geom_point(data = metadata %>% filter(Celltype %in% c("Sender", "Receiver")),
             aes(x = x, y = y, color = Celltype), size = 1.8, alpha = 0.9) +
  scale_color_manual(values = c("Sender" = "#990099", "Receiver" = "#0000FF")) +
  coord_fixed() +          
  labs(x = "x (\u00b5m)", y = "y (\u00b5m)", color = "Celltype",
       title = "Spatial distribution of Sender / Receiver cells") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),      # removes the grey gridlines
    panel.border = element_blank(),    # optional: drop the box border too, for a cleaner look
    axis.line = element_line(color = "grey40", linewidth = 0.4),
    plot.title = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))  # bigger, easier-to-see legend dots

ggsave(filename = plot_allneighbors_path, plot = p, width = 200, height = 150, units = "mm")

# diagnostic plot: illustrate the radii used by the CCC methods (via a single
# representative method, CellChat)
config = yaml::read_yaml("config.yaml")

vec_radius = config[["radius_param"]]$CellChat[[dataset]]

sender_data = metadata[metadata$Celltype == "Sender", ] %>% sample_n(50)

p = ggplot() +
  geom_point() +
  xlab("x") +
  ylab("y") +
  scale_color_manual(values = c("#990099", "#0000FF", "orange")) +
  scale_size_manual(values = c(2, 2, 1)) +
  theme_bw() +
  
  # Layer 1: largest radius
  geom_circle(data = sender_data,
              aes(x0 = x, y0 = y, r = vec_radius[4]),
              fill = "#DEEBF7", color = "#DEEBF7",
              alpha = 0.1, size = 0.5, linetype = "dashed", inherit.aes = FALSE) +
  
  # Layer 2
  geom_circle(data = sender_data,
              aes(x0 = x, y0 = y, r = vec_radius[3]),
              fill = "#9ECAE1", color = "#9ECAE1",
              alpha = 0.1, size = 0.5, linetype = "dashed", inherit.aes = FALSE) +
  
  # Layer 3
  geom_circle(data = sender_data,
              aes(x0 = x, y0 = y, r = vec_radius[2]),
              fill = "#4292C6", color = "#4292C6",
              alpha = 0.1, size = 0.5, linetype = "dashed", inherit.aes = FALSE) +
  
  # Layer 4: smallest radius
  geom_circle(data = sender_data,
              aes(x0 = x, y0 = y, r = vec_radius[1]),
              fill = "#084594", color = "#084594",
              alpha = 0.4, size = 0.5, linetype = "solid", inherit.aes = FALSE) +
  
  # Points on top
  geom_point(data = metadata,
             aes(x = x, y = y, color = Celltype, size = Celltype)) +
  
  scale_color_manual(values = c("#990099", "#0000FF", "orange")) +
  scale_size_manual(values = c(1.5, 1.5, 0.25))

ggsave(filename = plot_radius, plot = p, width = 200, height = 150, units = "mm")

# ==============================================================================
# STEP 6: estimate per-gene mean and dispersion (edgeR)
#
# Estimated using all neighbor pairs for each ligand (Sender) spot. These
# parameters drive the negative-binomial expression sampling used later in
# processing_semisimulation_nbsampling.R.
# ==============================================================================

mm = model.matrix(as.formula("~0 + Celltype"), metadata)

estimated_params = estimate_params_edgeR(counts = counts, metadata = metadata, mm = mm)

genemetadata = list(disp = data.frame(gene = rownames(estimated_params$dge), edgeR_dispersion = estimated_params$dge$tagwise.dispersion),
                    mean = estimated_params$means_perCT %>% as.data.frame,
                    edgeR_offsets = data.frame(cell_ID = colnames(counts), edgeR_offset = estimated_params$offset))

# ==============================================================================
# STEP 7: save outputs
# ==============================================================================

write.table(counts, processed_counts_path, sep = "\t")
saveRDS(genemetadata, genemetadata_path)
rownames(metadata) = NULL # remove rownames, otherwise the JSON file gets an extra column

write_json(list(metadata = metadata, neighbor_cells = neighbor_cells,
                average_percentageCells_expressingLR = round(average_percentageCells_expressingLR, 2), 
                average_distance_SenderReceiver = neighbors_info$distance_Sender_toClosest_Receiver %>% mean), cellmetadata_path)