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

# An useful error if the argument is missing
if (is.null(snakemake@input[["counts"]]) | is.null(snakemake@input[["metadata"]])){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

#############
### INPUT ###
#############
counts_path = snakemake@input[["counts"]]
metadata_path = snakemake@input[["metadata"]]

##############
### OUTPUT ###
##############
processed_counts_path = snakemake@output[["processed_counts"]]
genemetadata_path = snakemake@output[["genemetadata"]]
cellmetadata_path = snakemake@output[["cellmetadata"]]
plot_allneighbors_path = snakemake@output[["plot_allneighbors"]]
plot_radius = snakemake@output[["plot_radius"]]
plot_CT1CT2distance = snakemake@output[["plot_CT1CT2distance"]]
diagnostic_plots_path = snakemake@output[["diagnostic_plots"]]

##############
### PARAMS ###
##############
LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)
dataset = snakemake@params["dataset"] %>% unlist

###### Load data
counts = read.table(counts_path, row.names = 1)
metadata = read.table(metadata_path,row.names = 1)

# run QC pipeline without any kind of filtering
x = general_ST_QCpipeline_R(counts, to_filter = FALSE)

pdf(diagnostic_plots_path)
x[-1]
dev.off()

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
###### sample randomly CT1 and CT2 cells
metadata$Celltype = "Other"
metadata[sample(metadata%>% filter(Celltype == "Other") %>% pull(Cell_ID), 750), "Celltype"] = "CT1"
metadata[sample(metadata%>% filter(Celltype == "Other") %>% pull(Cell_ID), 750), "Celltype"] = "CT2"

###### Downsample datasets for celltypes besides CT1 and CT2
#set.seed(1)
#perc_dataset_to_remove = 0.5 
#cts = unique(metadata$Celltype) %>% setdiff(c("CT1","CT2"))
#for(celltype in cts)
#{
#  index_in_metadata = which(celltype == metadata$Celltype)
#  sampled_cells = sample(index_in_metadata, size = length(index_in_metadata) * perc_dataset_to_remove)
#  metadata %<>% filter(!row_number() %in% sampled_cells)
#  counts %<>% select(metadata$Cell_ID)
#}

# remove genes with 0 counts and keep genes expressed in at least 10 cells
#counts = counts[rowSums(counts) != 0 & rowSums(counts != 0) > 10,]

### STRATEGY 3 - manually select cells to create colocalization
# the removal of expression happens in the semisimukation script (search for STRATEGY 3)
### I use shiny code from script /CCC_benchmark_ST/shiny_forSelectingCells to manually select cells in space
#set.seed(1)
#if(dataset == "Visium_HD_HPC") {cells = c(read.csv("data/Visium_HD_HPC/DG.csv") %>% pull(Cell_ID) , read.csv("data/Visium_HD_HPC/CAx.csv") %>% pull(Cell_ID) )
#} else if(dataset == "MERFISH_mColon") {cells = read.csv("data/MERFISH_mColon/selected_cells.csv") %>% pull(Cell_ID)
#}
# split half cells to be CT1 and half CT2
#ct1cells = sample(cells, length(cells)*0.5)
#ct2cells = setdiff(cells, ct1cells)

#metadata[, "Celltype"] = "Other"
#metadata$Celltype[metadata$Cell_ID %in% ct1cells]= "CT1"
#metadata$Celltype[metadata$Cell_ID %in% ct2cells] = "CT2"

############################################################
#### Calculate percentage of cells expressing L/R genes  ###
############################################################
# We want to estimate how many cells actually express ligands (L) and receptor (R) genes
# estimate this for every gene and then average.
# this value will be the starting point for simulating percentage os cells/spots expressing LR genes

LRgenes = c(LR_database$ligand, LR_database$receptor) %>% 
  str_split("_") %>% 
  unlist %>% 
  unique()

average_percentageCells_expressingLR = apply(counts[which(rownames(counts) %in% LRgenes),] ,1,function(x) {
  sum(x > 0) / length(x)
}) %>% mean


######################################################
#### Select sender and neighboring receiver cells  ###
######################################################

# manually set values
# this values also correspond to config
# example config MERFISH_mColon: [40,75,100] 
#   - at 40 we should see anything as the minimum distance to receiver is 52
#   - at 75 we should pick all signal as all receivers are seen by methods
#   - at 100 we are checking how methods behave even tough they see same amount of interaction
if(dataset == "Visium_HD_HPC") {
  distance = 200
  distance_post_filtering_lower = 0 # 160
  distance_post_filtering_upper = 200
} else if(dataset == "MERFISH_mColon") {
  distance = 50
  distance_post_filtering_lower = 0 #42
  distance_post_filtering_upper = 75
}  else if(dataset == "CosMx_HFC") {
  distance = 200
  distance_post_filtering_lower = 0 # 150
  distance_post_filtering_upper = 250
} 

# there is a big variability for distances before and after sampling
# as we are trying to evaluate the distance parameter, we have to be more precise with sampling
# distance_threshold_dataset is needed because the function determines closest neighbor for each cell that must be above this value
# but if one calculates the closest neighbor after sampling, it is not exactly at distance_threshold_dataset as some cells are closer to others
# distance_post_filtering is important here as we will filter all connections of pairs of cells that are above/below this threshold
neighbors_info = find_neighboring_spots(spatial_coords = metadata %>% select(c("x","y")), 
                                        ligand_spots = metadata %>% filter(Celltype == "CT1") %>% select(Cell_ID) %>% unlist %>% unname, 
                                        receptor_spots = metadata %>% filter(Celltype == "CT2") %>% select(Cell_ID) %>% unlist %>% unname,
                                        remove_spots = TRUE,
                                        distance_threshold_dataset = distance,
                                        distance_post_filtering_lower = distance_post_filtering_lower,
                                        distance_post_filtering_upper = distance_post_filtering_upper)


table(neighbors_info$metadata$Celltype)

# remove from counts
counts = counts[,metadata$Cell_ID]

# Assuming your vectors are x and y
df = data.frame(d = neighbors_info$distance_CT1_toClosest_CT2)

# Create the main plot with fixed coordinates
p = ggplot(df, aes(x = seq_along(d), y = d)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(x = "Index of CT1", y = "Distance", title = "Distance to Closest CT2")


ggsave(filename = plot_CT1CT2distance, plot = p, width = 200, height = 150, units = "mm")


# ideally neighbors_info$average_distance_CT1_CT2_during_sampling and neighbors_info$average_distance_CT1_CT2_after_sampling should be very close to each other
# because for high FC of sender receiver cells, we want to test importance of l and average_distance_CT1_CT2_after_sampling is good indication for closes distance
# but for low FC, if average_distance_CT1_CT2_during_sampling is too high, then cells wont see each other at all and metrics wont reflect the reality


# Visium_HD_HPC (distance = 100) -> neighbors_info$average_distance_CT1_CT2_during_sampling -> 156 -> neighbors_info$average_distance_CT1_CT2_after_sampling -> 143
# MERFISH_mColon (distance = 50) -> neighbors_info$average_distance_CT1_CT2_during_sampling -> 52 -> neighbors_info$average_distance_CT1_CT2_after_sampling -> 43


########################
# In order to have same amount of cells across each strategy, sample cells according to var
#N_cells = 100

### STRATEGY 1- signal spread in space
# Update metadata according to find_neighboring_spots function
metadata = neighbors_info$metadata
neighbor_cells = neighbors_info$neighbors
# to standardize amount of cells that we have across STRATEGY 1 / 2 / 3 select only randomly N_cells cells
#neighbor_cells = neighbor_cells[,sample(neighbor_cells, N_cells) %>% names]

#metadata[, "Celltype"] = "Other"
#metadata[names(neighbor_cells), "Celltype"] = "CT1"
#metadata[unlist(neighbor_cells), "Celltype"] = "CT2"

### STRATEGY 2 - colocalization
### CHANGE CELLMETADATA TO ONLY CONTAIN CT1 CT2 IN SPECIFIC PLACE WHERE THEY ARE MORE ABUNDANT
#lst_out = select_HighDensityRegion( metadata = metadata, neighbor_cells = neighbor_cells, dataset = "Visium_HD_HPC", sample_cells = N_cells)

#neighbor_cells = lst_out$neighbor_cells
#metadata = lst_out$metadata

# simple plot
p = ggplot(metadata, aes(x = x, y = y,color = Celltype, size = Celltype)) +
  geom_point() +
  xlab("x") +
  ylab("y") +  
  scale_color_manual(values = c("#990099","#0000FF", "orange"))+
  scale_size_manual(values = c(2,2,1)) +
  theme_bw()



ggsave(filename = plot_allneighbors_path, plot = p, width = 200, height = 150, units = "mm")

## Create plot with radius to be used for methods
# as we try to standardize the amount of cells each method sees, here we select random method for illustration

config = yaml::read_yaml("config.yaml")

vec_radius = config[["l_param"]]$CellChat[[dataset]]

ct1_data = metadata[metadata$Celltype == "CT1", ] %>% sample_n(50)

p = ggplot() +
  theme_bw() +
  coord_fixed() + 
  
  # Layer 1: Largest
  geom_circle(data = ct1_data, 
              aes(x0 = x, y0 = y, r = vec_radius[4]), 
              fill = "#DEEBF7", color = "#DEEBF7", 
              alpha = 0.1, size = 0.5, linetype = "dashed", inherit.aes = FALSE) +
  
  # Layer 2
  geom_circle(data = ct1_data, 
              aes(x0 = x, y0 = y, r = vec_radius[3]), 
              fill = "#9ECAE1", color = "#9ECAE1", 
              alpha = 0.1, size = 0.5, linetype = "dashed", inherit.aes = FALSE) +
  
  # Layer 3
  geom_circle(data = ct1_data, 
              aes(x0 = x, y0 = y, r = vec_radius[2]), 
              fill = "#4292C6", color = "#4292C6", 
              alpha = 0.1, size = 0.5, linetype = "dashed", inherit.aes = FALSE) +
  
  # Layer 4: Smallest
  geom_circle(data = ct1_data, 
              aes(x0 = x, y0 = y, r = vec_radius[1]), 
              fill = "#084594", color = "#084594", 
              alpha = 0.4, size = 0.5, linetype = "solid", inherit.aes = FALSE) +
  
  # Points on top
  geom_point(data = metadata, 
             aes(x = x, y = y, color = Celltype, size = Celltype)) +
  
  scale_color_manual(values = c("#990099", "#0000FF", "orange")) +
  scale_size_manual(values = c(1.5, 1.5, 0.25))

#print("here")
#ggsave(filename = plot_radius, plot = p, width = 200, height = 150, units = "mm")

tiff(plot_radius, units="in", width=10, height=10, res=300)
p
dev.off()

######################################
# estimate mean and disp using edgeR #
######################################

# Parameters are estimated using all 6 neighbors for each ligand spot
mm= model.matrix(as.formula("~0 + Celltype") , metadata)

estimated_params = estimate_params_edgeR(counts = counts, metadata = metadata, mm = mm)

genemetadata = list( disp = data.frame(gene = rownames(estimated_params$dge) , edgeR_dispersion = estimated_params$dge$tagwise.dispersion) ,
                     mean = estimated_params$means_perCT %>% as.data.frame )


### save data
write.table(counts, processed_counts_path , sep = "\t")
saveRDS(genemetadata, genemetadata_path)
rownames(metadata) = NULL # remove rownames otherwise json file creates extra column
write_json(list(metadata = metadata, neighbor_cells = neighbor_cells, 
                average_percentageCells_expressingLR = round(average_percentageCells_expressingLR,2)),
           average_distance_CT1CT2 = neighbors_info$average_Practical_distance_CT1CT2, cellmetadata_path)
