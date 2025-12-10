'
Here I am using simple seurat wilcoxon statistical test to see how robustly it detects L-R interactions.
The goal here is not to perform QC, processing and clustering but rather immediately start statistical analysis assuming that I already know the clusters I am interested in
'

# Load package
suppressMessages({
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(jsonlite)
  library(magrittr)
  library(purrr)
  source("scripts/helper_functions.R")
})

# An useful error if the argument is missing
if (is.null(snakemake@input[["normalized_counts"]]) | is.null(snakemake@input[["cellmetadata_post_simulation"]]) 
    | is.null(snakemake@output[["significant_interactions"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read the argument
normalized_counts_path = snakemake@input[["normalized_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata_post_simulation"]]
significant_interactions_path = snakemake@output[["significant_interactions"]]

##############
### Params ###
##############
config = yaml::read_yaml("config.yaml")

dataset = snakemake@params["dataset"] %>% as.character

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)
#############
# load data #
#############

inflated_counts = read.csv(normalized_counts_path,sep="\t") %>% as.matrix
cellmetadata = read_json(path = cellmetadata_path)

'
inflated_counts = read.csv("output/Visium_HD_HPC_semiSimulation_NB//inflated_normalized_counts_FC_1_FC_nSenderCells_0.5_FC_nReceiverCells_1_indexLR_1.tsv",sep="\t") %>% as.matrix
cellmetadata = read_json(path = "output/Visium_HD_HPC_semiSimulation_NB//simulated_cellmetadata_FC_1_FC_nSenderCells_0.5_FC_nReceiverCells_1_indexLR_1.json")
LR_database = read.table("data/LR_database.tsv", row.names = 1)
'

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

# Create seurat object
# The seurat boject will contain only normalized counts
SO_obj = CreateSeuratObject(counts = inflated_counts, meta.data = cellmetadata$metadata)
SO_obj@assays$RNA$data = inflated_counts

# Seurat only does pairwise comparisons without subunit info. This means that cases like L-R1-R2 are not possible
# To bypass this, I split these cases into L-R1 and L-R2. The analysis results are averaged

# Filter LR_database and Seurat object to only contain common ligand receptor genes
# Transform LR database from L-R1-R2 to L-R1 | L-R2
lst = apply(LR_database, 1, as.list)

LR_database = lapply(lst, function(x) {
  grid = x[c("ligand","receptor")] %>%
    str_split(.,"_") %>%
    expand.grid() %>%
    mutate(ligand_receptor = x$ligand_receptor) %>%
    dplyr::rename(ligand = Var1,
           receptor = Var2,)
  
  all_present = (grid %>% select("ligand","receptor") %>% unlist ) %in% rownames(SO_obj) %>%
    all
  
  if(all_present) {return(grid)}
}) %>% Reduce(rbind.data.frame,.) 

# remove duplicated entries as some interactions are the same (ex: WNT2B_FZD7_LRP6 and WNT2B_FZD4_LRP6)
n = LR_database %>% 
  select(ligand,receptor) %>% 
  apply(.,1, function(x) {str_c(x,collapse = "_")}) %>% 
  duplicated %>% 
  which

LR_database = LR_database[-n,]

# filter SO
SO_obj = SO_obj[c(LR_database$ligand,LR_database$receptor) %>% unique %>% as.character,]

# Wilcoxon test
# some genes wont be present in output because of min.cells.feature/min.cells.group parameters
# Perform differential expression to find ligands that are overexpressed in CT1
Idents(SO_obj) = SO_obj$Celltype
markers = FindMarkers(SO_obj, slot = "data" , ident.1 = "CT1", ident.2 = "CT2", test.use = "wilcox")
markers_ligands = markers %>% filter(avg_log2FC > 0)
# Perform differential expression to find receptors that are overexpressed in CT2
markers = FindMarkers(SO_obj, slot = "data" , ident.1 = "CT2", ident.2 = "CT1", test.use = "wilcox")
markers_receptors = markers %>% filter(avg_log2FC > 0)

# Average the results for cases like L-R1-R2 (aggregating L-R1 and L-R2)
seurat_LRR_averaged_out = list()
for(entry in unique(LR_database$ligand_receptor))
{
  name = str_split(entry,"_") %>% 
    unlist
  
  ligand = name[1] # 1st entry is alwasy ligand
  if(length(name) == 2) {
    receptor = name[2]
    
    # Average over the results
    markers_ligands_subset = markers_ligands %>% filter(rownames(.) %in% ligand)
    markers_receptors_subset = markers_receptors %>% filter(rownames(.) %in% receptor)
    
    # if at least one of the genes in "name" is not present in the output, then the aggregated result doesnt exist
    if(nrow(markers_ligands_subset) != 1 | nrow(markers_receptors_subset) != 1) {next}
    
    tmp_markers = rbind(markers_ligands_subset, markers_receptors_subset)
    averaged_results = tmp_markers %>%
      summarise(
        across(where(is.character), ~ first(.)),
        across(where(is.factor), ~ first(.)), 
        across(where(is.numeric), ~ mean(., na.rm = TRUE)))
    
    seurat_LRR_averaged_out[[entry]] = averaged_results
    
  } else if (length(name) == 3){
    receptor = c(name[2], name[3])
    
    # Average over the results
    markers_ligands_subset = markers_ligands %>% filter(rownames(.) %in% ligand)
    markers_receptors_subset = markers_receptors %>% filter(rownames(.) %in% receptor)
    
    # if at least one of the genes in "name" is not present in the output, then the aggregated result doesnt exist
    if(nrow(markers_ligands_subset) != 1 | nrow(markers_receptors_subset) != 2) {next}
    
    tmp_markers = rbind(markers_ligands_subset, markers_receptors_subset)
    averaged_results = tmp_markers %>%
      summarise(
        across(where(is.character), ~ first(.)),
        across(where(is.factor), ~ first(.)), 
        across(where(is.numeric), ~ mean(., na.rm = TRUE)))
    
    seurat_LRR_averaged_out[[entry]] = averaged_results
    
  } else {message("seurat_wilcoxon analysis complex with more than 2 receptors") ; break}
    
 
} 
seurat_LRR_averaged_out = do.call(rbind.data.frame,seurat_LRR_averaged_out) %>%
  mutate(ligand_receptor = rownames(.),
         significant = p_val_adj < 0.05,
         statistics = p_val_adj) %>%
  select(p_val, p_val_adj , ligand_receptor,significant,statistics) %>% # drop the log2fc and pct columns as they dont reflect reality
  arrange(statistics) # sort from lower to higher pvalues


# save data
write.table(seurat_LRR_averaged_out ,significant_interactions_path)