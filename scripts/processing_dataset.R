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
diagnostic_plots_path = snakemake@output[["diagnostic_plots"]]

##############
### PARAMS ###
##############
LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)

###### Load data
counts = read.table(counts_path, row.names = 1)
metadata = read.table(metadata_path,row.names = 1)

# run QC pipeline without any kind of filtering
x = general_ST_QCpipeline_R(counts, to_filter = FALSE)

pdf(diagnostic_plots_path)
x[-1]
dev.off()

'
counts = read.table("data/Visium_HD_HPC/counts_Visium_HD_HPC.tsv", row.names = 1)
metadata = read.table("data/Visium_HD_HPC/metadata_Visium_HD_HPC.tsv",row.names = 1)
'

###### Downsample datasets 2x
set.seed(1)
perc_dataset_to_remove = 0.5 # variable specify percentage of data to randomly remove to speed up workflow
for(celltype in unique(metadata$Celltype))
{
  index_in_metadata = which(celltype == metadata$Celltype)
  sampled_cells = sample(index_in_metadata, size = length(index_in_metadata) * perc_dataset_to_remove) # remove 50% of sample
  metadata %<>% filter(!row_number() %in% sampled_cells)
  counts %<>% select(metadata$Cell_ID)
}

# remove genes with 0 counts and keep genes expressed in at least 10 cells
counts = counts[rowSums(counts) != 0 & rowSums(counts != 0) > 10,]

### STRATEGY 3 - manually select cells to create colocalization
# the removal of expression happens in the semisimukation script (search for STRATEGY 3)
### I use shiny code from script /CCC_benchmark_ST/shiny_forSelectingCells to manually select cells in space
#set.seed(1)
cells = c(read.csv("data/Visium_HD_HPC/DG.csv") %>% pull(Cell_ID) , read.csv("data/Visium_HD_HPC/CAx.csv") %>% pull(Cell_ID) )
# split half cells to be CT1 and half CT2
ct1cells = sample(cells, length(cells)*0.5)
ct2cells = setdiff(cells, ct1cells)

metadata[, "Celltype"] = "Other"
metadata[ct1cells, "Celltype"] = "CT1"
metadata[ct2cells, "Celltype"] = "CT2"

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
neighbors_info = find_neighboring_spots(spatial_coords = metadata %>% select(c("x","y")), 
                                        ligand_spots = metadata %>% filter(Celltype == "CT1") %>% select(Cell_ID) %>% unlist %>% unname, 
                                        receptor_spots = metadata %>% filter(Celltype == "CT2") %>% select(Cell_ID) %>% unlist %>% unname,
                                        remove_spots = FALSE)


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
                average_percentageCells_expressingLR = round(average_percentageCells_expressingLR,2)), cellmetadata_path)
