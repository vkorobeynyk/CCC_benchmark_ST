suppressMessages({
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(edgeR)
  library(ggpubr)
  library(sf)
  library(jsonlite)
  source("scripts/helper_functions.R")
})

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
FC_nSenderCells = snakemake@wildcards[["FC_nSenderCells"]] %>% as.double
FC_nReceiverCells = snakemake@wildcards[["FC_nReceiverCells"]] %>% as.double
LR_database_path = snakemake@params[["LR_database"]]
max_N_neighbors = snakemake@params[["max_N_neighbors"]] %>% as.integer
#################
### Load data ###
#################

counts = read.table(processed_counts_path)
genemetadata = readRDS(genemetadata_path)
means_perCT = genemetadata$mean
cellmetadata = read_json(path = cellmetadata_path)

'
counts = read.table("data/processed/Visium_HD_HPC/processed_counts_Visium_HD_HPC.tsv")
genemetadata = readRDS("data/processed/Visium_HD_HPC/genemetadata_Visium_HD_HPC.RDS")
means_perCT = genemetadata$mean
cellmetadata = read_json("data/processed/Visium_HD_HPC/cellmetadata_Visium_HD_HPC.json")
LRdb = read.table("data/LR_database.tsv", header = T)
'

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

######################
### Important info ###
######################
# cellmetadata contains all celltypes and all cells
# In the end of the script, I am saving the cellmetadata as the metadata that cell-cell communications methods will use

# check if cell names of counts and metadata correspond and are in the same order
stopifnot(colnames(counts) == cellmetadata$metadata$Cell_ID)

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
# We change the amount of receiver/sender cells that express receptor 
fraction_cells_expressingR = FC_nReceiverCells * cellmetadata$average_percentageCells_expressingLR
fraction_cells_expressingL = FC_nSenderCells * cellmetadata$average_percentageCells_expressingLR

# Semi simulation
# Subunits are also simulated because the LR database nomenclature is L_R1_R2 etc
semi_simulation_out = semi_simulate(counts = counts, simulated_interactions_lst = simulated_interactions_lst , genemetadata = genemetadata, 
                                    metadata = cellmetadata$metadata ,
                                    fraction_cells_expressingR = fraction_cells_expressingR ,
                                    fraction_cells_expressingL = fraction_cells_expressingL,
                                    FC = FC,
                                    df_neighbors = cellmetadata$neighbor_cells)

# simple plot to show cells to which we added signal
cellmetadata$metadata$Celltype_updated = cellmetadata$metadata$Celltype
cellmetadata$metadata$Celltype_updated[cellmetadata$metadata$Cell_ID %in% semi_simulation_out$cell_info$CT1_CT2$cells_signalAdded$CT1] = "CT1_signalAdded"
cellmetadata$metadata$Celltype_updated[cellmetadata$metadata$Cell_ID %in% semi_simulation_out$cell_info$CT1_CT2$cells_signalAdded$CT2] = "CT2_signalAdded"

# simple plotfind_neighboring_spots
p = ggplot(cellmetadata$metadata, aes(x = x, y = y,color = Celltype_updated, size = Celltype_updated)) +
  geom_point() +
  xlab("x") +
  ylab("y")+  
  scale_color_manual(values = c("#FFCCFF" ,"#990099" ,"#CCCCFF" ,"#0000FF" ,"#FFCC99")) +
  scale_size_manual(values = c(1.5,3,1.5,3,0.5))

ggsave(filename = plot_neighbors_path, plot = p, width = 200, height = 150, units = "mm")

################
# save results #
################
write.table(semi_simulation_out$counts_inflated, inflated_counts_path , sep = "\t")
saveRDS(simulated_interactions_lst,simulated_interactions_path)
#saveRDS(FC_after_semisimulation,FC_after_simulation_path)

# replace metadata and neighbors information with newly computed neighbors according to *find_neighboring_spots* function
# as those cells expression were modified
cellmetadata$neighbor_cells = cellmetadata$neighbor_cells
write_json(cellmetadata, simulated_cellmetadata_path)
#write_json(list(cellmetadata$metadata, neighbor_cells = cellmetadata$neighbor_cells, cell_coordinates = metadata[,c("x","y")]), simulated_cellmetadata_path)

### Cellphonedb strictly requires specific metadata file
## Create here
mt = cellmetadata$metadata %>% select(c("Cell_ID", "Celltype")) %>% setNames(c("barcode_sample","cell_type"))
rownames(mt) = mt$barcode_sample
write.table(mt, metadata_cpdbv5_path , sep = "\t")