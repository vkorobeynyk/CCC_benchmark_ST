# ============================================================================
# Runs CellChat on the semi-simulated spatial dataset and extracts the
# significance/probability of the simulated Sender -> Receiver ligand-receptor
# interaction.
# ============================================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(CellChat)
  library(jsonlite)
  library(magrittr)
  library(purrr)
  source("scripts/helper_functions.R")
})

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

radius_index = as.integer(snakemake@params["radius_index"]) + 1
dataset = snakemake@params["dataset"] %>% as.character
radius = config[["radius_param"]][["CellChat"]][[dataset]][radius_index] %>% unlist

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)

# ==============================================================================
# STEP 1: load data
# ==============================================================================

inflated_counts = read.csv(normalized_counts_path, sep = "\t") %>% as.matrix
cellmetadata = read_json(path = cellmetadata_path)

# ---------------------------------------------------------------------------
# Manual interactive testing snippet (not executed): quickly load a single
# run's files in an R session, e.g. for debugging outside Snakemake.
# ---------------------------------------------------------------------------
'
inflated_counts = read.csv("output/Visium_HD_HPC_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_0.5_PCE_Receiver_6_indexLR_2.tsv",sep="\t") %>% as.matrix
cellmetadata = read_json(path = "output/Visium_HD_HPC_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_0.5_PCE_Receiver_6_indexLR_2.json")
LR_database = read.table("data/LR_database.tsv", row.names = 1)
'

# transform the JSON list into individual dataframes
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

# ==============================================================================
# STEP 2: figure out which cells the method actually "sees"
#
# Build the pairwise euclidean distance matrix, then for every Sender cell
# that received inflated signal, find every other cell within `radius` of it.
# ==============================================================================

coord = cellmetadata$metadata %>%
  select(x, y) %>%
  mutate_at(vars(x, y), as.numeric)

distance_mat = apply(coord, 1, function(pt)
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
  ggtitle("Cellchat euclidean filtering pink -> Sender | orange -> cells seen by method | purple -> Sender_signalAdded | blue -> Receiver_signalAdded") +
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
    ratio_Receiver_seen_byMethod = amount_Receiver_seen_byMethod/ amount_Sender_signalAdded_cells * 100
  } else if (amount_Sender_signalAdded_cells == amount_Receiver_signalAdded_cells) {
    ratio_Receiver_seen_byMethod = amount_Receiver_seen_byMethod / amount_Receiver_signalAdded_cells * 100
  }
} else {
  ratio_Receiver_seen_byMethod = 0
}

# average number of cells that each Sender cell has "seen" by the method
average_cells_perSender_seen_byMethod = length(all_cells_seen_byMethod) / amount_Sender_signalAdded_cells

# ==============================================================================
# STEP 4: build the CellChat object and restrict its LR database to LR_database.tsv
# ==============================================================================

# CellChat, if needed, transforms coordinates from pixels to micrometers using 2
# parameters: ratio and tol. Ratio is the pixel-to-micrometer conversion rate;
# since all our data is already in micrometers, ratio = 1. Tol is usually half of
# the cell/spot size, and matters when comparing center-to-center distance
# against the "interaction range" parameter; since we just want every method to
# use the same spatial distance, tol isn't important here and is set to 1.
# All datasets were rescaled to micrometers.
spatial.factors = data.frame(ratio = 1, tol = 1)

cellchat = createCellChat(object = inflated_counts, meta = data.frame(Celltype = cellmetadata$metadata$Celltype,
                                                                      samples = "sample1" %>% as.factor,
                                                                      row.names = cellmetadata$metadata$Cell_ID),
                          group.by = "Celltype", datatype = "spatial", coordinates = as.matrix(coord), spatial.factors = spatial.factors)

CellChatDB = CellChatDB.human

# ---------------------------------------------------------------------------
# To understand how CellChatDB works, we looked at the example (from
# LR_database.tsv) of:
#   ligand (L) -> WNT4
#   receptor (R) -> LRP6_FZD8
# CellChatDB has the interaction as WNT4_FZD8_LRP6, but LR_database has
# WNT4_LRP6_FZD8. Here we switch the receptor order (we cannot switch the
# order in CellChatDB directly, otherwise CellChat would no longer recognize
# the interaction as valid).
# ---------------------------------------------------------------------------

fix_gene_order = function(df) {
  sapply(df, function(x)
  {
    parts = strsplit(x, "_")[[1]]
    if(length(parts) == 2) {str_c(parts[2] , "_" , parts[1])
    } else {x}
  })
}

LR_database$receptor = fix_gene_order(LR_database$receptor)
LR_database$ligand = fix_gene_order(LR_database$ligand)
LR_database$interaction_name  = str_c(LR_database$ligand, "_", LR_database$receptor)

# CellChatDB is a list with 4 entries:
#  - interaction: subset to the same interactions as LR_database.tsv
#  - geneInfo: per-gene information, no need to filter
#  - complex: no need to filter, CellChat only fetches what it needs
#  - cofactor: we aren't simulating cofactors, so we simply blank them out
CellChatDB$interaction = CellChatDB$interaction[which(CellChatDB$interaction$interaction_name %in% LR_database$interaction_name), ]

stopifnot(nrow(CellChatDB$interaction) == nrow(LR_database))

x = CellChatDB$cofactor %>% apply(., 2, function(x) {return(rep("", length(x)))})
rownames(x) = rownames(CellChatDB$cofactor)
CellChatDB$cofactor = x %>% as.data.frame() # CellChat requires a dataframe here

cellchat@DB = CellChatDB # same dim() as LR_database

cellchat = subsetData(cellchat) # necessary even when using the whole database

# ==============================================================================
# STEP 5: run CellChat and extract the significance/probability of the
# Sender -> Receiver interaction
# ==============================================================================

# Wilcoxon test to filter features down (only the p-value threshold is used here)
cellchat = identifyOverExpressedGenes(cellchat, min.cells = 0, thresh.fc = 0, thresh.p = 0.05)
cellchat = identifyOverExpressedInteractions(cellchat)

# if the Wilcoxon test didn't find any significant LR pair to test spatially
if (nrow(cellchat@LR$LRsig) == 0)
{
  write.table(data.frame(ligand_receptor = NA, significant = FALSE, statistics = 0,
                         ratio_Receiver_seen_byMethod = ratio_Receiver_seen_byMethod,
                         average_cells_perSender_seen_byMethod = average_cells_perSender_seen_byMethod), significant_interactions_path,
              sep = "\t")
} else {
  
  # Since we're testing how many cells should express each gene, trim = 0.001 means
  # only 0.1% of cells need to express the gene to be counted
  # When comparing communication across different CellChat objects, the same scale factor should be used
  # I tested different scale.distance (0.1,0.5,1) and the results didnt change. What changes is the probability but pvalue always stays the same
  cellchat = computeCommunProb(cellchat, type = "truncatedMean", trim = 0.001,
                               distance.use = TRUE, interaction.range = radius, scale.distance = 1.5 / (distance_mat[distance_mat > 0] %>% min), # 1.5 midpoint of [1,2]
                               contact.dependent = FALSE, contact.range = NULL, nboot = 100)
  
  df.net = subsetCommunication(cellchat, thresh = 1) # p-value threshold for determining a significant interaction
  
  df.net %<>% filter(source == "Sender" & target == "Receiver") %>% mutate(ligand_receptor = gsub("—", "_", interaction_name),
                                                                           significant = pval < 0.05, # only returns significant interactions
                                                                           statistics = prob) %>% dplyr::arrange(desc(prob))
  
  # save data
  if (nrow(df.net) != 0)
  {
    write.table(data.frame(ligand_receptor = df.net$interaction_name, significant = df.net$significant, statistics = df.net$prob,
                           ratio_Receiver_seen_byMethod = ratio_Receiver_seen_byMethod,
                           average_cells_perSender_seen_byMethod = average_cells_perSender_seen_byMethod), significant_interactions_path,
                sep = "\t")
  } else {
    write.table(data.frame(ligand_receptor = NA, significant = FALSE, statistics = 0,
                           ratio_Receiver_seen_byMethod = ratio_Receiver_seen_byMethod,
                           average_cells_perSender_seen_byMethod = average_cells_perSender_seen_byMethod), significant_interactions_path,
                sep = "\t")
  }
}
