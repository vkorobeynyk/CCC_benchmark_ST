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

###### Load data
counts = read.table(counts_path, row.names = 1)
metadata = read.table(metadata_path,row.names = 1)

# run QC pipeline without any kind of filtering
x = general_ST_QCpipeline_R(counts, to_filter = FALSE)

pdf(diagnostic_plots_path)
x[-1]
dev.off()

'
counts = read.table("data/STARmap_plus_HPC/counts_STARmap_plus_HPC.tsv", row.names = 1)
metadata = read.table("data/STARmap_plus_HPC/metadata_STARmap_plus_HPC.tsv",row.names = 1)
'

###### Downsample datasets by 2x
set.seed(1)
for(celltype in unique(metadata$Celltype))
{
  index_in_metadata = which(celltype == metadata$Celltype)
  sampled_cells = sample(index_in_metadata, size = length(index_in_metadata) * 1/2)
  metadata %<>% filter(!row_number() %in% sampled_cells)
  counts %<>% select(metadata$Cell_ID)
}


# remove genes with 0 counts and keep genes expressed in at least 10 cells
counts = counts[rowSums(counts) != 0 & rowSums(counts != 0) > 10,]

######################################################
#### Select sender and neighboring receiver cells  ###

#n_neighbors = 6 -> this is used to estimate mean and dispersion. As there is no perfect number for number of neighbors, we chose one that gives enough interactions 
# that can be simulated

neighbors_info = find_neighboring_spots(spatial_coords = metadata %>% select(c("x","y")), 
                                        n_neighbors = 6,  
                                        ligand_spots = metadata %>% filter(Celltype == "CT1") %>% select(Cell_ID) %>% unlist %>% unname, 
                                        receptor_spots = metadata %>% filter(Celltype == "CT2") %>% select(Cell_ID) %>% unlist %>% unname,
                                        remove_spots = TRUE)

# simple plot
p = ggplot(neighbors_info$metadata, aes(x = x, y = y,color = Celltype, size = Celltype)) +
  geom_point() +
  xlab("x") +
  ylab("y") +  
  scale_color_manual(values = c("#008000","orange", "black"))+
  scale_size_manual(values = c(2,2,0.75))

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
write_json(list(metadata = metadata, neighbor_cells = neighbors_info$neighbors), cellmetadata_path)
