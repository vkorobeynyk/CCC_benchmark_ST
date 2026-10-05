# ============================================================================
# Runs NICHES on the semi-simulated spatial dataset and extracts the
# significance of the simulated Sender -> Receiver ligand-receptor interaction.
# ============================================================================

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(NICHES)
  library(jsonlite)
  library(magrittr)
  library(purrr)
})
source("scripts/helper_functions.R")

# Fail early with a clear error if a required input/output is missing
if (is.null(snakemake@input[["normalized_counts"]]) | is.null(snakemake@input[["cellmetadata_post_simulation"]])
    | is.null(snakemake@output[["plot_neighbors"]]) | is.null(snakemake@output[["significant_interactions"]])) {
  stop("Argument_name needs to be specified, but is missing.n", call. = FALSE)
}

# ==============================================================================
# STEP 0: Snakemake I/O and params
# ==============================================================================

# input files
normalized_counts_path = snakemake@input[["normalized_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata_post_simulation"]]

# output files
significant_interactions_path = snakemake@output[["significant_interactions"]]
plot_neighbors_path = snakemake@output[["plot_neighbors"]]

# params
config = yaml::read_yaml("config.yaml")

radius_index = (as.integer(snakemake@params["radius_index"]) + 1)
dataset = snakemake@params["dataset"] %>% as.character
radius = config[["radius_param"]][["NICHES"]][[dataset]][radius_index] %>% unlist

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)

# ==============================================================================
# STEP 1: load data
# ==============================================================================

inflated_counts = read.csv(normalized_counts_path, sep = "\t") %>% as.matrix
cellmetadata = read_json(path = cellmetadata_path)

# transform the JSON list into individual dataframes
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

SO = CreateSeuratObject(counts = inflated_counts, assay = "logcounts", meta.data = cellmetadata$metadata)
SO@assays$logcounts$data = SO@assays$logcounts$counts # NICHES uses the "data" slot by default, with no argument to change this

# ==============================================================================
# STEP 2: figure out which cells the method actually "sees"
# ==============================================================================

coord = cellmetadata$metadata %>%
  select(x, y) %>%
  mutate_at(vars(x, y), as.numeric)

distance_mat <- apply(coord, 1, function(pt)
  (sqrt(abs(pt["x"] - coord$x)^2 + abs(pt["y"] - coord$y)^2))
)

Sender = cellmetadata$metadata %>% filter(Celltype == "Sender") %>% pull(Cell_ID)
Receiver = cellmetadata$metadata %>% filter(Celltype == "Receiver") %>% pull(Cell_ID)
Sender_signalAdded = cellmetadata$metadata %>% filter(Celltype_updated == "Sender_signalAdded") %>% pull(Cell_ID)
Receiver_signalAdded = cellmetadata$metadata %>% filter(Celltype_updated == "Receiver_signalAdded") %>% pull(Cell_ID)

# for each Sender_signalAdded cell, find all cells within `radius` of it
vec = map(Sender_signalAdded, function(cell_OI) {
  
  within_radius = distance_mat[, cell_OI == colnames(distance_mat)] < radius
  return(which(within_radius))
  
}) %>% as.list

all_cells_seen_byMethod = colnames(distance_mat)[unlist(vec) %>% unique]
all_cells_seen_byMethod = all_cells_seen_byMethod[!all_cells_seen_byMethod %in% Sender_signalAdded] # remove Sender cells themselves

# diagnostic plot: which cells fall within the method's neighborhood radius
plt = ggplot(coord, aes(x = x, y = y)) +
  geom_point(size = 0.1) +
  geom_point(data = coord[all_cells_seen_byMethod, ], aes(x = x, y = y), colour = "orange", size = 2) +
  geom_point(data = coord[Sender, ], aes(x = x, y = y), colour = "#FFCCFF", size = 2) +
  geom_point(data = coord[Sender_signalAdded, ], aes(x = x, y = y), colour = "#990099", size = 3) +
  geom_point(data = coord[Receiver_signalAdded, ], aes(x = x, y = y), colour = "#0000FF", size = 3) +
  ggtitle("NICHES euclidean filtering pink -> Sender | orange -> cells seen by method | purple -> Sender_signalAdded | blue -> Receiver_signalAdded") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme_bw()

ggsave(plot_neighbors_path, plt, device = "png", width = 30, height = 25, units = "cm")

# ==============================================================================
# STEP 3: compute neighborhood-coverage metrics
#
# How many of the Receiver_signalAdded cells actually fall within the
# method's neighborhood definition, and how many cells on average does each
# Sender cell "see"?
# ==============================================================================

amount_Sender_signalAdded_cells = cellmetadata$metadata %>% filter(Celltype_updated == "Sender_signalAdded") %>% nrow
amount_Receiver_signalAdded_cells = cellmetadata$metadata %>% filter(Celltype_updated == "Receiver_signalAdded") %>% nrow
amount_Receiver_seen_byMethod = which(all_cells_seen_byMethod %in% (cellmetadata$metadata %>% filter(Celltype_updated == "Receiver_signalAdded") %>% pull(Cell_ID))) %>% length

# PCE_Sender > PCE_Receiver -> how many Receivers are seen by the method?
# PCE_Sender < PCE_Receiver -> do all Senders see at least 1 Receiver?
# PCE_Sender == PCE_Receiver -> are all Receivers seen by both Sender and the method?
if (amount_Receiver_seen_byMethod != 0) {
  if (amount_Sender_signalAdded_cells > amount_Receiver_signalAdded_cells) {
    ratio_Receiver_seen_byMethod = amount_Receiver_seen_byMethod / amount_Receiver_signalAdded_cells * 100
  } else if (amount_Sender_signalAdded_cells < amount_Receiver_signalAdded_cells) {
    ratio_Receiver_seen_byMethod =  amount_Receiver_seen_byMethod/ amount_Sender_signalAdded_cells * 100
  } else if (amount_Sender_signalAdded_cells == amount_Receiver_signalAdded_cells) {
    ratio_Receiver_seen_byMethod = amount_Receiver_seen_byMethod / amount_Receiver_signalAdded_cells * 100
  }
} else {
  ratio_Receiver_seen_byMethod = 0
}

# average number of cells that each Sender cell has "seen" by the method
average_cells_perSender_seen_byMethod = length(all_cells_seen_byMethod) / amount_Sender_signalAdded_cells

# ==============================================================================
# STEP 4: run NICHES and extract markers of the Sender -> Receiver niche
# ==============================================================================

# add metadata variables to the Seurat object
SO@meta.data %<>% mutate(Celltype = cellmetadata$metadata$Celltype,
                         x = coord$x,
                         y = coord$y)

NICHES_output = RunNICHES(object = SO,
                          LR.database = "custom",
                          custom_LR_database = LR_database,
                          assay = "logcounts",
                          position.x = "x",
                          position.y = "y",
                          k = NULL,
                          rad.set = radius, # spatial radius to account for
                          cell_types = "Celltype",
                          min.cells.per.ident = 0,
                          min.cells.per.gene = NULL,
                          meta.data.to.map = c('Celltype'),
                          CellToCell = F, CellToSystem = F, SystemToCell = F,
                          CellToCellSpatial = T, CellToNeighborhood = F, NeighborhoodToCell = F)

niche_CtN = NICHES_output[['CellToCellSpatial']]

# Wilcoxon test for markers of the Sender -> Receiver niche
markers_CtN = FindAllMarkers(niche_CtN, test.use = "wilcox") %>% filter(cluster == "Sender—Receiver")

markers_CtN %<>% mutate(ligand_receptor = gsub("—|-", "_", gene),
                        significant = markers_CtN$p_val_adj < 0.05,
                        statistics = p_val_adj) %>%
  arrange(p_val_adj)

# ==============================================================================
# STEP 5: save results
# ==============================================================================

if (nrow(markers_CtN) != 0)
{
  markers_CtN$ratio_Receiver_seen_byMethod = ratio_Receiver_seen_byMethod
  markers_CtN$average_cells_perSender_seen_byMethod = average_cells_perSender_seen_byMethod
  write.table(markers_CtN, significant_interactions_path, sep = "\t")
} else {
  write.table(data.frame(ligand_receptor = NA, significant = FALSE, statistics = 0,
                         ratio_Receiver_seen_byMethod = ratio_Receiver_seen_byMethod,
                         average_cells_perSender_seen_byMethod = average_cells_perSender_seen_byMethod), significant_interactions_path, sep = "\t")
}