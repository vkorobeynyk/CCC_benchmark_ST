# Load package
library(mistyR)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(jsonlite)
source("scripts/helper_functions.R")

# An useful error if the argument is missing
if (is.null(snakemake@input[["processed_counts"]]) | is.null(snakemake@input[["cellmetadata"]]) | is.null(snakemake@output[["cellmetadata2"]])){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}
# Read the argument
processed_counts_path = snakemake@input[["processed_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata"]]
cellmetadata2_path = snakemake@output[["cellmetadata2"]]
plot_neighbors_path = snakemake@output[["plot_neighbors"]]

#############
# load data #
#############
processed_counts = read.table(processed_counts_path) %>% as.matrix
cellmetadata = read_json(cellmetadata_path)

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)

rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

'
processed_counts = read.table("data/processed/CosMx_HFC/processed_counts_CosMx_HFC.tsv") %>% as.matrix
cellmetadata = read_json("data/processed/CosMx_HFC/cellmetadata3_CosMx_HFC.json")
'

coord = cellmetadata$metadata %>% select(x,y)

amount_neighbors_seen = 1
l = 0.1
while(any(amount_neighbors_seen != 6))
{
  amount_neighbors_seen = vector()
  
  l = l+0.1
  for(cell_OI in colnames(cellmetadata$neighbor_cells))
  {
    n = which(rownames(coord) == cell_OI) # get index of cell
    cell_OI_coord = coord %>% dplyr::slice(n)
    neighbor_cells = cellmetadata$neighbor_cells[cell_OI] %>% unlist
    cells_all = rownames(coord)
    
    # calculate euclidean distance from Cell_OI to all other cells
    eucl_dist_vec = sqrt((rep(cell_OI_coord[,"x"],length(cells_all)) - coord[cells_all,"x"])^2) + sqrt((rep(cell_OI_coord[,"y"],length(cells_all)) - coord[cells_all,"y"])^2)
    
    # calculate weights based on formula from mistyR docs for gaussian
    coord$eucl_dist = eucl_dist_vec
    w = exp(-(eucl_dist_vec^2 / l^2)) 
    coord$w = w
    
    # check that all neighbor cells that we added signal to are seen by the method
    amount_neighbors_seen = append(amount_neighbors_seen, sum(coord[neighbor_cells,"w"] > 0.1))
  }
}

# save figure showing sender , receiver and all cells with positive weights of last cell
plt = ggplot(coord, aes(x = x ,y = y, color = w)) + 
  geom_point() +
  geom_point(data=coord[neighbor_cells,] , aes(x=x, y=y), colour="red", size=2) +
  ggtitle("log transformed weights based on formula exp(-(eucl_dist_vec^2 / l^2))")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

ggsave(plot_neighbors_path, plt, device = "png", width = 30, height = 25, units = "cm")

# save the spatial weight parameter
cellmetadata$spatialWeight_mistyR = l

# save file
write_json(cellmetadata, cellmetadata2_path)
