'
NICHES allows to look at cell-cell comunication at the single cell level or at the niche level.
At single cell level -> NICHES builds a matrix of rows (L-R pair) and columns (pairs of all cells). THen this matrix can be used as
normal gene expression matrix, that can be clustered, to find communication patterns between cells.
At niche level (CelltoNeighboorhood) -> NICHES builds matrix of rows (L-R pair) and columns (sending cell - neighborhood according to radius).
One can do analysis like at single cell level and perform differential expression analysis to find what communication patterns Celltype1
sends that other celltypes do not send.
'

# Load package
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(NICHES)
library(rjson)
library(magrittr)
library(purrr)
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

radius_index = (as.integer(snakemake@params["radius_index"]) +1)
dataset = snakemake@params["dataset"] %>% as.character
radius = config[["l_param"]][["NICHES"]][[dataset]][radius_index] %>% unlist

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)
#############
# load data #
#############

inflated_counts = read.csv(normalized_counts_path,sep="\t") %>% as.matrix
cellmetadata = fromJSON(file = cellmetadata_path)

'
inflated_counts = read.csv("output/CosMx_HFC_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_2_indexLR_26.tsv",sep="\t") %>% as.matrix
cellmetadata = fromJSON(file = "output/CosMx_HFC_semiSimulation_NB/simulated_cellmetadata_1_n_neigbors_2_indexLR_26.json")
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

# generate a list where each index name is sender cell and it contains all cells within the radius seen by the method
all_sender_cells = colnames(cellmetadata$neighbor_cells)
all_receiver_cells = unlist(cellmetadata$neighbor_cells)
vec = map(all_sender_cells, function(cell_OI) {
  simulated_neighbors = cellmetadata$neighbor_cells[,cell_OI]
  
  within_radius = distance_mat[,grep(cell_OI, colnames(distance_mat))] < radius
  return(within_radius)
  
}) %>% as.data.frame()

rownames(vec) = colnames(distance_mat)
all_cells_seen_byMethod = vec[rowSums(vec)>0,] %>% rownames

# plot
plt = ggplot(coord, aes(x = x ,y = y)) + 
  geom_point(size = 0.5) +
  geom_point(data=coord[all_cells_seen_byMethod,] , aes(x=x, y=y), colour="orange", size=2) +
  geom_point(data=coord[all_receiver_cells,] , aes(x=x, y=y), colour="#008000", size=2) +
  geom_point(data=coord[all_sender_cells,] , aes(x=x, y=y), colour="blue", size=2) +
  ggtitle("NICHES euclidean radius filtering | orange -> cells within radius | blue -> sender cells | green -> receiver cells")+
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
markers_CtN = FindAllMarkers(niche_CtN,min.pct = 0,test.use = "wilcox") %>% filter(cluster == "CT1")

markers_CtN %<>% mutate(ligand_receptor = gsub("—","_",rownames(markers_CtN)),
                        significant = markers_CtN$p_val < 0.05,
                        statistics = p_val) %>% 
  arrange(p_val)
# save data
write.table(markers_CtN ,significant_interactions_path)