'
NICHES allows to look at cell-cell comunication at the single cell level or at the niche level.
At single cell level -> NICHES builds a matrix of rows (L-R pair) and columns (pairs of all cells). THen this matrix can be used as
normal gene expression matrix, that can be clustered, to find communication patterns between cells.
At niche level (CelltoNeighboorhood) -> NICHES builds matrix of rows (L-R pair) and columns (sending cell - neighborhood according to radius).
One can do analysis like at single cell level and perform differential expression analysis to find what communication patterns Celltype1
sends that other celltypes do not send.
'
suppressMessages({
  # Load package
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

# An useful error if the argument is missing
if (is.null(snakemake@input[["normalized_counts"]]) | is.null(snakemake@input[["cellmetadata_post_simulation"]]) 
    | is.null(snakemake@output[["plot_neighbors"]]) | is.null(snakemake@output[["significant_interactions"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read the argument
normalized_counts_path = snakemake@input[["normalized_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata_post_simulation"]]
significant_interactions_path = snakemake@output[["significant_interactions"]]
plot_neighbors_path = snakemake@output[["plot_neighbors"]]

##############
### Params ###
##############
config = yaml::read_yaml("config.yaml")

l_index = (as.integer(snakemake@params["l_index"]) +1)
dataset = snakemake@params["dataset"] %>% as.character
radius = config[["l_param"]][["NICHES"]][[dataset]][l_index] %>% unlist

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)
#############
# load data #
#############

inflated_counts = read.csv(normalized_counts_path,sep="\t") %>% as.matrix
cellmetadata = read_json(path = cellmetadata_path)

'
inflated_counts = read.csv("output/Visium_HD_HPC_semiSimulation_NB//inflated_normalized_counts_FC_1_FC_nSenderCells_7_FC_nReceiverCells_0.5_indexLR_1.tsv",sep="\t") %>% as.matrix
cellmetadata = read_json(path = "output/Visium_HD_HPC_semiSimulation_NB//simulated_cellmetadata_FC_1_FC_nSenderCells_7_FC_nReceiverCells_0.5_indexLR_1.json")
'

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

# Create Seurat object
SO = CreateSeuratObject(counts = inflated_counts, assay = "logcounts", meta.data = cellmetadata$metadata)
SO@assays$logcounts$data = SO@assays$logcounts$counts # NICHES uses "data" slot by befaul and doesnt have argument to change

########################################################
# Compute all sender and receiver cells in the dataset #
########################################################

# CODE ADAPTED FROM NICHES GITHUB TO CALCULATE THE RADIUS FOR PLOTTING REASONS
# Compute the euclidean distance matrix
coord = cellmetadata$metadata %>% 
  select(x,y) %>% 
  mutate_at(vars(x,y) , as.numeric)

distance_mat <- apply(coord, 1, function(pt)
  (sqrt(abs(pt["x"] - coord$x)^2 + abs(pt["y"] - coord$y)^2))
)
rownames(distance_mat) = colnames(distance_mat)
# generate a list where each index name is sender cell and it contains all cells within the radius seen by the method
CT1 = colnames(cellmetadata$neighbor_cells)
CT2 = unlist(cellmetadata$neighbor_cells)[1]
CT1_signalAdded = cellmetadata$metadata$Cell_ID[cellmetadata$metadata$Celltype_updated == "CT1_signalAdded"]
CT2_signalAdded = cellmetadata$metadata$Cell_ID[cellmetadata$metadata$Celltype_updated == "CT2_signalAdded"]
vec = map(CT1, function(cell_OI) {
  within_radius = distance_mat[,which(cell_OI == colnames(distance_mat))] < radius
  return(which(within_radius))
  
}) %>% as.list

all_cells_seen_byMethod = lapply(vec, names) %>% unlist

# plot
plt = ggplot(coord, aes(x = x ,y = y)) + 
  geom_point(size = 0.1) +
  geom_point(data=coord[all_cells_seen_byMethod,] , aes(x=x, y=y), colour="orange", size=2) +
  geom_point(data=coord[CT1_signalAdded,] , aes(x=x, y=y), colour="#990099", size=3) +
  geom_point(data=coord[CT2_signalAdded,] , aes(x=x, y=y), colour="#0000FF", size=3) +
  ggtitle("NICHES euclidean radius filtering orange -> cells seen by method | purple -> CT1_signalAdded | blue -> CT2_signalAdded")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

ggsave(plot_neighbors_path, plt, device = "png", width = 30, height = 25, units = "cm")

# Add metadata variables to SO
SO@meta.data %<>% mutate(Celltype = cellmetadata$metadata$Celltype,
                         x = coord$x,
                         y = coord$y)

##############
# Run method #
##############
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
                          CellToCell = F,CellToSystem = F,SystemToCell = F,
                          CellToCellSpatial = F,CellToNeighborhood = T,NeighborhoodToCell = F)

niche_CtN = NICHES_output[['CellToNeighborhood']]

# CellToNeighborhood analysis (here we have signal averaging over k neighbors)
# Find cell-cell communication between CT1 and neighborhood comparing to CT2-Neighborhood and Other-Neighborhood
markers_CtN = FindAllMarkers(niche_CtN, test.use = "wilcox") %>% filter(cluster == "CT1")

markers_CtN %<>% mutate(ligand_receptor = gsub("—","_",gene),
                        significant = markers_CtN$p_val_adj < 0.05,
                        statistics = p_val_adj) %>% 
  arrange(p_val_adj)
# save data
write.table(markers_CtN ,significant_interactions_path)