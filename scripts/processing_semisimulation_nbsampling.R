library(dplyr)
library(stringr)
library(liana)
library(magrittr)
library(edgeR)
library(ggpubr)
library(sf)
library(jsonlite)
source("scripts/helper_functions.R")

# An useful error if the argument is missing
if (is.null(snakemake@input[["processed_counts"]]) | is.null(snakemake@input[["genemetadata"]]) | is.null(snakemake@params[["LR_database"]]) | 
    is.null(snakemake@input[["cellmetadata"]])){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}
#############
### INPUT ###
#############
processed_counts_path = snakemake@input[["processed_counts"]]
genemetadata_path = snakemake@input[["genemetadata"]]
cellmetadata_path = snakemake@input[["cellmetadata"]]

##############
### OUTPUT ###
##############
inflated_counts_path = snakemake@output[["inflated_counts"]]
simulated_cellmetadata_path = snakemake@output[["simulated_cellmetadata"]]
simulated_interactions_path = snakemake@output[["simulated_interactions"]]
FC_after_simulation_path = snakemake@output[["FC_after_simulation"]]
metadata_cpdbv5_path = snakemake@output[["metadata_cpdbv5"]]
plot_neighbors_path = snakemake@output[["plot_neighbors"]]

##############
### PARAMS ###
##############
indexLR_toSample = snakemake@wildcards[["indexLR_toSample"]] %>% as.integer
FC = snakemake@wildcards[["FC"]] %>% as.double
n_neighbors = snakemake@wildcards[["n_neighbors"]] %>% as.integer
FC_nCells_expressingLR = snakemake@wildcards[["FC_nCells_expressingLR"]] %>% as.double
LR_database_path = snakemake@params[["LR_database"]]

#################
### Load data ###
#################

counts = read.table(processed_counts_path)
genemetadata = readRDS(genemetadata_path)
means_perCT = genemetadata$mean
cellmetadata = read_json(cellmetadata_path)

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

######################
### Important info ###
######################
# there might be confusion with metadata
# cellmetadata contains all celltypes and all cells
# neighbors_metadata is the output of find_neighboring_spots and contains all spots (CT1 and CT2) whose counts will be inflated
# neighbors_metadata is used to indicate cells where the signal will be added and to estimate parameters (using edgeR) for calculation of FC after semi-simulation
# In the end of the script, I am saving the cellmetadata as the metadata that cell-cell communications methods will use



'
counts = read.table("data/processed/MERFISH_mColon//processed_counts_MERFISH_mColon.tsv")
genemetadata = readRDS("/home/vkorob/Documents/git/CCC_benchmark_ST/data/processed/MERFISH_mColon/genemetadata_MERFISH_mColon.RDS")
means_perCT = genemetadata$mean
cellmetadata = read_json("data/processed/MERFISH_mColon/cellmetadata_MERFISH_mColon.json")
LRdb = read.table("data/LR_database.tsv", header = T)
'

#######################################
# Select n_neighbors to be simulated  #
#######################################
#### Select neighboring cells / spots according to n_neighbors parameter
# This is not the final amounht of cells to which the signal will be added to
# The final amount will be sampled from all neighbors selected above * FC_nCells_expressingLR (probability of a cell expressing Receptor)
neighbors_info = cellmetadata$neighbor_cells[1:n_neighbors,]

# Add Celltype information to metadata file
neighbors_metadata = cellmetadata$metadata %>% mutate(Celltype = ifelse(Cell_ID %in% colnames(neighbors_info) , "CT1", "Other"))
neighbors_metadata$Celltype[neighbors_metadata$Cell_ID %in% unlist(unname(neighbors_info)) & neighbors_metadata$Celltype != "CT1"] = "CT2"

# check if cell names of counts and metadata correspond and are in the same order
stopifnot(colnames(counts) == neighbors_metadata$Cell_ID)

message(paste("Shape of counts object:" , str_flatten(dim(counts) , " ")))

########################################################
# remove genes with mean == 0 in celltypes to simulate #
########################################################

gene_index = which(means_perCT[,"CT1"] == 0 | means_perCT[,"CT2"] == 0)
if(length(gene_index) > 0) {
  means_perCT = means_perCT[-gene_index,]
}

#############
# Load LRdb #
#############

# Load database
LRdb = read.table(LR_database_path, header = T)

# iterate over every row and filter genes that are not in count data and have 0 mean per CT
n = LRdb %>%
  apply(., 1, function(row) {
    str_split(row,"_") %>%
      lapply(., function(x) {(x %in% rownames(counts)) & (x %in% means_perCT$gene_names)}) %>%
      unlist %>%
      all
})
LRdb = LRdb[n,]

message(paste("After filtering LR database," , nrow(LRdb) , "LR pairs show expression in at least 10 cells"))

###########################################
# Select genes to artificially add counts #
###########################################
simulated_interactions_lst = list()

# Pre-sample LR pairs to be used in the semi simulation
# Also pre-sample the subunit genes
comb_CT = "CT1_CT2"

# select LR pair to inflate expression according to index indexLR_toSample
LR_sample = LRdb[indexLR_toSample,]
simulated_interactions_lst[[comb_CT]]$ligand = LR_sample$ligand
simulated_interactions_lst[[comb_CT]]$receptor = LR_sample$receptor

###########################
# Inflate gene expression #
###########################
# calculate the percentage of cells to which add signal to
# For the case of senders, we keep the amount of cells to which add expression to the same through the entire benchmark
# We change the amount of receiver cells that express receptor 
# This is merely to decrease computational cost as adding another parameter would increase the computational time dramatically
N_cells_expressingR = FC_nCells_expressingLR * 100 * cellmetadata$average_percentageCells_expressingLR
N_cells_expressingL = 100 * cellmetadata$average_percentageCells_expressingLR

# Semi simulation
# Subunits are also simulated because the LR database nomenclature is L_R1_R2 etc
semi_simulation_out = semi_simulate(counts = counts, simulated_interactions_lst = simulated_interactions_lst , genemetadata = genemetadata, 
                                    metadata = neighbors_metadata ,
                                    N_cells_expressingR = N_cells_expressingR ,
                                    N_cells_expressingL = N_cells_expressingL,
                                    FC = FC,
                                    df_neighbors = neighbors_info)

# simple plot to show cells to which we added signal
neighbors_metadata$Celltype[neighbors_metadata$Cell_ID %in% semi_simulation_out$cell_info$CT1_CT2$cells_signalAdded$CT1] = "CT1_signalAdded"
neighbors_metadata$Celltype[neighbors_metadata$Cell_ID %in% semi_simulation_out$cell_info$CT1_CT2$cells_signalAdded$CT2] = "CT2_signalAdded"

# simple plotfind_neighboring_spots
p = ggplot(neighbors_metadata, aes(x = x, y = y,color = Celltype, size = Celltype)) +
  geom_point() +
  xlab("x") +
  ylab("y")+  
  scale_color_manual(values = c("#FFCCFF" ,"#990099" ,"#CCCCFF" ,"#0000FF" ,"#FFCC99")) +
  scale_size_manual(values = c(1.5,3,1.5,3,0.5))

ggsave(filename = plot_neighbors_path, plot = p, width = 200, height = 150, units = "mm")

#################################
# calculate FC after simulation #
#################################
# estimate the mean only for cells that had their expression increased
# as the cells that we plan to add signal to were determined with *find_neighboring_spots* function. THere are no cells to which we didnt add signal

# get simulated genes
names_ofLgenes = simulated_interactions_lst$CT1_CT2$ligand %>% str_split(., "_")  %>% unlist
names_ofRgenes = simulated_interactions_lst$CT1_CT2$receptor %>% str_split(., "_")   %>% unlist

FC_after_semisimulation = list()
for(CT in names(semi_simulation_out$cell_info[[1]]$cells_signalAdded))
{
  # select the proper set of genes (either ligand for sender CT or receptor for receiving CT)
  if(CT == names(semi_simulation_out$cell_info[[1]]$cells_signalAdded)[1]) {simulated_genes = names_ofLgenes} else if(CT == names(semi_simulation_out$cell_info[[1]]$cells_signalAdded)[2]) {simulated_genes = names_ofRgenes}
  
  stopifnot(neighbors_metadata$Cell_ID == colnames(semi_simulation_out$counts_inflated)) # just for security
  
  # create model matrix
  mm= model.matrix(as.formula("~0 + Celltype") , neighbors_metadata)
  
  # Estimate parameters
  estimate_params_out = estimate_params_edgeR(counts = semi_simulation_out$counts_inflated, metadata = neighbors_metadata,mm = mm)
  
  FC_after_semisimulation[[CT]] = estimate_params_out$means_perCT %>% arrange(., gene_names) %>% filter(.,gene_names %in% simulated_genes) %>% .[CT] / 
    means_perCT %>% arrange(., gene_names) %>% filter(.,gene_names %in% simulated_genes) %>% .[CT]
}


################
# save results #
################
write.table(semi_simulation_out$counts_inflated, inflated_counts_path , sep = "\t")
saveRDS(simulated_interactions_lst,simulated_interactions_path)
saveRDS(FC_after_semisimulation,FC_after_simulation_path)

# replace metadata and neighbors information with newly computed neighbors according to *find_neighboring_spots* function
# as those cells expression were modified
cellmetadata$metadata = neighbors_metadata
cellmetadata$neighbor_cells = neighbors_info
write_json(cellmetadata, simulated_cellmetadata_path)
#write_json(list(cellmetadata$metadata, neighbor_cells = cellmetadata$neighbor_cells, cell_coordinates = metadata[,c("x","y")]), simulated_cellmetadata_path)

### Cellphonedb strictly requires specific metadata file
## Create here
mt = cellmetadata$metadata %>% select(c("Cell_ID", "Celltype")) %>% setNames(c("barcode_sample","cell_type"))
rownames(mt) = mt$barcode_sample
write.table(mt, metadata_cpdbv5_path , sep = "\t")