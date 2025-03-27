# Load package
library(mistyR)
library(magrittr)
library(stringr)
library(jsonlite)
library(ggplot2)
library(dplyr)
source("scripts/helper_functions.R")

# An useful error if the argument is missing
if (is.null(snakemake@input[["normalized_counts"]]) | is.null(snakemake@input[["cellmetadata_post_simulation"]]) | is.null(snakemake@output[["significant_interactions"]]) |
    is.null(snakemake@params[["LR_database"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}
# Read the argument
normalized_counts_path = snakemake@input[["normalized_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata_post_simulation"]]
LR_database_path = snakemake@params[["LR_database"]]
significant_interactions_path = snakemake@output[["significant_interactions"]]
results_folder = snakemake@output[["results_folder"]]

#############
# load data #
#############
normalized_counts = read.table(normalized_counts_path) %>% as.matrix
cellmetadata = read_json(cellmetadata_path)
LRdb = read.table(LR_database_path, header = T)

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

##########
# PARAMS #
##########
l = cellmetadata$spatialWeight_mistyR %>% unlist
  
'
normalized_counts = read.table("output/STARmap_plus_HPC_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_4.tsv") %>% as.matrix
cellmetadata = read_json("output/STARmap_plus_HPC_semiSimulation_NB/simulated_cellmetadata_1_n_neigbors_4.json")
LRdb = read.table("data/LR_database.tsv", header = T)

'
##############
# Run method #
##############
# mistyR build a framework composed of different views. Each view will test individual genes
# In order to speed up the computation time, I select only genes that are ligands and receptors as the metrics are calculated only based on those

# Select L and R 
n = which(rownames(normalized_counts) %in% LRdb$ligand)
n1 = which(rownames(normalized_counts) %in% LRdb$receptor) 
normalized_counts = normalized_counts[union(n,n1),] 

normalized_counts %<>% t %>% as.data.frame() # mistyR needs genes as columns

# misty doesnt like genes with "-" -> replace here and after misty calculations
colnames(normalized_counts) = gsub("-","_",colnames(normalized_counts))

misty.intra = create_initial_view(normalized_counts)
summary(misty.intra)
summary(misty.intra$intraview)

coord = cellmetadata$metadata %>% 
  select(x,y) %>%
  mutate(across(c(x,y),as.numeric))

#The juxtaview represent a local spatial view and captures the expression of all markers available in the intraview within the immediate neighborhood of each cell. 
#The paraview captures the expression of all markers available in the intraview in the boarder tissue structure where the importance of the influence is proportional to the inverse of the distance between two cells
misty.views = misty.intra %>% add_paraview(., positions = coord, l = l, family = "gaussian")
#summary(misty.views)

misty.views %>% run_misty(.,results.folder=results_folder, cached = FALSE)

misty.results = collect_results(results_folder)

method_out = misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>% 
  summarize(mean.p = mean(value)) %>%
  arrange(mean.p) %>%
  as.data.frame

######################
# Create output list #
######################

# filter all genes that were found significant
method_out %<>% filter(mean.p < 0.05)

# transform back the "_" to "-"
method_out$target = gsub("_","-",method_out$target)

# remove rows where Ligand and Receptor are the same
# the line below remove cases where ligand and receptor are same gene. Makes sense as usually ligand cant be receptor (and vice-versa), and people doing analysis know that
method_out_filtered = method_out$target %>% expand.grid(.,.) %>% filter(Var1 != Var2) %>% rename(ligand = Var1, receptor = Var2)
method_out_filtered$ligand_receptor = str_c(method_out_filtered$ligand , "_", method_out_filtered$receptor)

# add pvalue column based on mean value of original pvalues -> this may be used for future ranking of interactions
pvalues_ligand = method_out[method_out_filtered$ligand,"mean.p"]
pvalues_receptors = method_out[method_out_filtered$receptor,"mean.p"]
method_out_filtered$pval = rowMeans(cbind(pvalues_ligand , pvalues_receptors))

method_out_filtered %<>% arrange(pval)

# save file
saveRDS(method_out_filtered , significant_interactions_path)
