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
inflated_counts = read.csv("output/MERFISH_mColon_semiSimulation_NB/inflated_normalized_counts_FC_1_FC_nSenderCells_4_FC_nReceiverCells_9_indexLR_1.tsv",sep="\t") %>% as.matrix
cellmetadata = read_json(path = "output/MERFISH_mColon_semiSimulation_NB/simulated_cellmetadata_FC_1_FC_nSenderCells_4_FC_nReceiverCells_9_indexLR_1.json")
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
lst = split(LR_database, seq_len(nrow(LR_database)))

LR_database = lapply(lst, function(x) {
  grid = x[c("ligand","receptor")] %>%
    str_split(.,"_") %>%
    expand.grid() %>%
    mutate(ligand_receptor = x$ligand_receptor) 
  
  #colnames(grid)[1:2] = c("ligand","receptor")
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
markers = FindAllMarkers(SO_obj, slot = "data", test.use = "wilcox")

# Initialize as empty dataframes with correct column names if no markers are found
if (nrow(markers) == 0) {
  markers_ligands = markers
  markers_receptors = markers
} else {
  markers_ligands = markers %>% filter(cluster == "CT1") %>% filter(avg_log2FC > 0)
  markers_receptors = markers %>% filter(cluster == "CT2") %>% filter(avg_log2FC > 0)
}

seurat_LRR_averaged_out = list()

# Process L-R Database
if (nrow(markers_ligands) > 0 & nrow(markers_receptors) > 0) {
  
  for(entry in unique(LR_database$ligand_receptor)) {
    name = unlist(str_split(entry, "_"))
    ligand = name[1]
    
    # Handle cases for 1 or 2 receptors (L_R or L_R1_R2)
    if(length(name) %in% c(2, 3)) {
      receptors = name[2:length(name)]
      
      # Subset markers for the specific pair
      subset_L = markers_ligands %>% filter(gene == ligand)
      subset_R = markers_receptors %>% filter(gene %in% receptors)
      
      # Check if ALL components are found
      # (Ligand must be 1, Receptors must match the count in the database string)
      if(nrow(subset_L) == 1 && nrow(subset_R) == (length(name) - 1)) {
        
        tmp_markers = rbind(subset_L, subset_R)
        
        averaged_results = tmp_markers %>%
          dplyr::summarise(
            across(where(is.character), ~ data.table::first(.)),
            across(where(is.factor), ~ data.table::first(.)), 
            across(where(is.numeric), ~ mean(., na.rm = TRUE))
          )
        
        seurat_LRR_averaged_out[[entry]] = averaged_results
      }
    } else {
      message(paste("Skipping complex pair:", entry))
    }
  }
}

if(length(seurat_LRR_averaged_out) > 0) {
  seurat_LRR_averaged_out = do.call(rbind.data.frame, seurat_LRR_averaged_out) %>%
    mutate(
      ligand_receptor = rownames(.),
      significant = p_val_adj < 0.05,
      statistics = p_val_adj
    ) %>%
    select(p_val, p_val_adj, ligand_receptor, significant, statistics) %>%
    arrange(statistics)
} else {
  # Default output if no L-R pairs were significantly co-expressed
  seurat_LRR_averaged_out = data.frame(
    p_val = 1, 
    p_val_adj = 1, 
    ligand_receptor = "none_detected", 
    significant = FALSE, 
    statistics = 1
  )
}

# Metadata columns
seurat_LRR_averaged_out$ratio_CT2_seen_byMethod = 100
seurat_LRR_averaged_out$average_cells_perCT1_seen_byMethod = FALSE

# save data
write.table(seurat_LRR_averaged_out ,significant_interactions_path)