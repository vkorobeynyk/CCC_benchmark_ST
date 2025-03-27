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
LR_database_path = snakemake@params[["LR_database"]]
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
nLR_per_CTCTcomb = snakemake@params[["nLR_per_CTCTcomb"]]
FC = as.double(snakemake@wildcards[["FC"]])
n_neighbors = as.integer(snakemake@wildcards[["n_neighbors"]])

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
counts = read.table("data/processed/STARmap_plus_HPC/processed_counts_STARmap_plus_HPC.tsv")
genemetadata = readRDS("/home/vkorob/Documents/git/CCC_benchmark_ST/data/processed/STARmap_plus_HPC/genemetadata_STARmap_plus_HPC.RDS")
means_perCT = genemetadata$mean
cellmetadata = read_json("data/processed/STARmap_plus_HPC/cellmetadata4_STARmap_plus_HPC.json")
LRdb = read.table("data/LR_database.tsv", header = T)
'

################################################
# Select appropriate amount of receptor spots  #
################################################
#### Select neighboring cells / spots acording to n_neighbors parameter
neighbors_info = cellmetadata$neighbor_cells[1:n_neighbors,]

# Add Celltype information to metadata file
neighbors_metadata = cellmetadata$metadata %>% mutate(Celltype = ifelse(Cell_ID %in% colnames(neighbors_info) , "CT1", "Other"))
neighbors_metadata$Celltype[neighbors_metadata$Cell_ID %in% unlist(unname(neighbors_info)) & neighbors_metadata$Celltype != "CT1"] = "CT2"

# simple plotfind_neighboring_spots
p = ggplot(neighbors_metadata, aes(x = x, y = y,color = Celltype, size = Celltype)) +
  geom_point() +
  xlab("x") +
  ylab("y")+  
  scale_color_manual(values = c("#0072B2","#D55E00", "#41DE11")) +
  scale_size_manual(values = c(2,2,0.75))

ggsave(filename = plot_neighbors_path, plot = p, width = 200, height = 150, units = "mm")

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

# Load OmniPath database
LRdb = read.table(LR_database_path, header = T)

LRdb = LRdb[(LRdb$ligand %in% rownames(counts)),]
LRdb = LRdb[(LRdb$receptor %in% rownames(counts)),]
# filter LR because there are duplicated pairs (only 1)
LRdb$L_R = str_c(LRdb$ligand,"_",LRdb$receptor)
LRdb = LRdb[!duplicated(LRdb$L_R),]

# filter LRdb to only contain L/R that have mean != 0
LRdb %<>% filter(ligand %in% means_perCT$gene_names)
LRdb %<>% filter(receptor %in% means_perCT$gene_names)

message(paste("After filtering LR database," , nrow(LRdb) , "LR pairs show expression in at least 10 cells"))

###########################
# Inflate gene expression #
###########################
set.seed(1)
simulated_interactions_lst = list()

tmp_LRdb = LRdb

# Pre-sample LR pairs to be used in the semi simulation
# Also pre-sample the subunit genes
# This is written in case we want to simulate multiple interactions and not only CT1_CT2
for(comb_CT in "CT1_CT2")
{
  # select randomly LR pair to inflate expression
  index_toSample = sample(seq(1,nrow(tmp_LRdb)) , size = nLR_per_CTCTcomb, replace = F)
  LR_sample = tmp_LRdb[index_toSample,]
  
  # Add genes belonging to a complex to be also inflated (based on CellPhoneDB and CellChatDB DB)
  # CellChat
  LRdb_cellchat = select_resource(c('CellChatDB'))[[1]]
  LRdb_cellchat = rbind(filter(LRdb_cellchat, grepl("COMPLEX", source)) , filter(LRdb_cellchat, grepl("COMPLEX", target)))
  
  # CellPhoneDB
  LRdb_cpdb = select_resource(c('CellPhoneDB'))[[1]]
  LRdb_cpdb = rbind(filter(LRdb_cpdb, grepl("COMPLEX", source)) , filter(LRdb_cpdb, grepl("COMPLEX", target)))
  
  LRdb_forSubunits_joined = rbind(LRdb_cellchat , LRdb_cpdb)
  LRdb_forSubunits_joined = LRdb_forSubunits_joined[!str_c(LRdb_forSubunits_joined$source, LRdb_forSubunits_joined$target) %>% duplicated,] # remove duplicated entries
  
  # remove the sampled LR pairs 
  tmp_LRdb = tmp_LRdb[-index_toSample,]
  
  # adding L and R to the LR_list to track inflated genes without removing duplicates
  simulated_interactions_lst[[comb_CT]] = LR_sample$L_R
  
  ############## find subunits for ligands
  for(gene in LR_sample$ligand)
  {
    tmp_df = filter(LRdb_forSubunits_joined , grepl(gene, source_genesymbol) & grepl("COMPLEX", source))
    
    # in case this gene has no subunits, skip
    if(nrow(tmp_df) == 0) {next}
    
    subunits = c(tmp_df$source_genesymbol %>% str_split("_") %>% lapply("[",1) , tmp_df$source_genesymbol %>% str_split("_") %>% lapply("[",2)) %>% unlist %>% unique()
    subunits = subunits[!grepl(gene, subunits)] # remove original gene
    
    # dont select subunits that are not in means_perCT (probably were filtered because estimated mean = 0 or the gene doesnt exist in the data)
    subunits %<>% .[subunits %in% means_perCT$gene_names]
    
    simulated_interactions_lst[[comb_CT]]  %<>% append(. , subunits[which(subunits %in% rownames(counts))] %>% str_c(., "_subunit")) # remove empty strings and add subunit . Also here we filter subunits that are not present in count data
  }
  
  ############## find subunits for receptors
  for(gene in LR_sample$receptor)
  {
    tmp_df = filter(LRdb_forSubunits_joined , grepl(gene, target_genesymbol) & grepl("COMPLEX", target))
    
    # in case this gene has no subunits, skip
    if(nrow(tmp_df) == 0) {next}
    
    subunits = c(tmp_df$target_genesymbol %>% str_split("_") %>% lapply("[",1) , tmp_df$target_genesymbol %>% str_split("_") %>% lapply("[",2)) %>% unlist %>% unique()
    subunits = subunits[!grepl(gene, subunits)] # remove original gene
    
    # dont select subunits that are not in means_perCT (probably were filtered because estimated mean = 0 or the gene doesnt exist in the data)
    subunits %<>% .[subunits %in% means_perCT$gene_names]
    
    simulated_interactions_lst[[comb_CT]]  %<>% append(. , subunits[which(subunits %in% rownames(counts))] %>% str_c("subunit_" , .)) # remove empty strings and add subunit . Also here we filter subunits that are not present in count data
  }
  
  # remove duplicates in subunits info
  simulated_interactions_lst[[comb_CT]] = simulated_interactions_lst[[comb_CT]][!duplicated(simulated_interactions_lst[[comb_CT]])]
}

# Semi simulation
# It may happen that a subunit of a gene has a mean parameter that is below the 1Q threshold I use. Now I use the original means for the subunits. In theory i would have 
# to check the parameters of every subunit and they are not in line, I would remove them
semi_simulation_out = semi_simulate(counts = counts, simulated_interactions_lst = simulated_interactions_lst , genemetadata = genemetadata, 
                                    metadata = neighbors_metadata , combination_CT = "CT1_CT2", FC = FC, n_neighbors = n_neighbors, df_neighbors = neighbors_info)

#################################
# calculate FC after simulation #
#################################
# estimate the mean only for cells that had their expression increased
# as the cells that we plan to add signal to were determined with *find_neighboring_spots* function. THere are no cells to which we didnt add signal

# get simulated genes
names_ofLgenes = str_split(simulated_interactions_lst[[1]], "_")  %>% lapply("[[",1) %>% unlist %>% setdiff(.,"subunit") %>% unique
names_ofRgenes = str_split(simulated_interactions_lst[[1]], "_")  %>% lapply("[[",2) %>% unlist %>% setdiff(.,"subunit") %>% unique

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