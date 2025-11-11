'
Giotto is an analysis toolbox for spatial transcriptomic data. The way it computes LR interaction is by first creating
a spatial graph of neighbor cells based on knn, distance cut-off or delanuay triangulation.
For each LR pair it calculates score S for every pair of cells. Latter is based on weighted average expression.
To assess significance it creates a random null distribution by shuffling cell locations and calculates pvalue based on frequency
of S from null distribution being greater or smaller to original distribution.
A final differential activity score is calculated by multiplying the log2 fold change with the adjusted p values
'

# Load package
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(Giotto)
  library(jsonlite)
  library(magrittr)
  library(purrr)
  source("scripts/helper_functions.R")
})

# An useful error if the argument is missing
if (is.null(snakemake@input[["normalized_counts"]]) | is.null(snakemake@input[["cellmetadata_post_simulation"]]) 
    | is.null(snakemake@output[["plot_neighbors"]]) | is.null(snakemake@output[["significant_interactions"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read the argument
normalized_counts_path = snakemake@input[["normalized_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata_post_simulation"]]
significant_interactions_path = snakemake@output[["significant_interactions"]]
plot_neighbors_path = snakemake@output[["plot_neighbors"]]

##############
### Params ###
##############
config = yaml::read_yaml("config.yaml")

l_index = (as.integer(snakemake@params["l_index"]) +1)
dataset = snakemake@params["dataset"] %>% as.character
radius = config[["l_param"]][["giotto"]][[dataset]][l_index] %>% unlist

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)
#############
# load data #
#############

inflated_counts = read.csv(normalized_counts_path,sep="\t") %>% as.matrix
cellmetadata = read_json(path = cellmetadata_path)

'
inflated_counts = read.csv("output/Visium_HD_HPC_semiSimulation_NB//inflated_normalized_counts_FC_1_FC_nSenderCells_0.5_FC_nReceiverCells_0.5_indexLR_1.tsv",sep="\t") %>% as.matrix
cellmetadata = read_json(path = "output/Visium_HD_HPC_semiSimulation_NB//simulated_cellmetadata_FC_1_FC_nSenderCells_0.5_FC_nReceiverCells_0.5_indexLR_1.json")
'

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID
cellmetadata$metadata %<>%
  dplyr::rename(cell_ID = Cell_ID) # giotto wants column name cell_ID

# Create Giotto object
giotto_obj = createGiottoObject(expression = inflated_counts,
                                   spatial_locs = cellmetadata$metadata %>% select(x,y),
                                 cell_metadata = cellmetadata$metadata)


# Giotto only do pairwise comparisons without subunit info. This means that cases like L-R1-R2 are not possible
# To bypass this, I split these cases into L-R1 and L-R2. The analysis results are averaged

# Filter LR_database and Giotto object to only contain common ligand receptor genes
# as giotto uses only expression of the latter genes, this doesnt influence the analysis
# Transform LR database from L-R1-R2 to L-R1 | L-R2
lst = apply(LR_database, 1, as.list)

LR_database = lapply(lst, function(x) {
  grid = x[c("ligand","receptor")] %>%
    str_split(.,"_") %>%
    expand.grid() %>%
    mutate(ligand_receptor = x$ligand_receptor) %>%
    dplyr::rename(ligand = Var1,
           receptor = Var2,)
  
  all_present = (grid %>% select("ligand","receptor") %>% unlist ) %in% rownames(giotto_obj) %>%
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

# filter giotto object
giotto_obj = giotto_obj[c(LR_database$ligand,LR_database$receptor) %>% unique %>% as.character,]

# Create spatial network
giotto_obj = createSpatialNetwork(gobject = giotto_obj,
                                    maximum_distance_delaunay = radius, 
                                    method = 'Delaunay', 
                                    delaunay_method = 'delaunayn_geometry')



# Giotto requires normalization to be performed
# Since we are already importing normalized data, replace in the giotto object
giotto_obj = suppressWarnings(normalizeGiotto(giotto_obj)) # suppress warning as it complains about library size
giotto_obj@expression$cell$rna$normalized = giotto_obj@expression$cell$rna$raw

# visualize the spatial network
plt = spatPlot(gobject = giotto_obj, show_network = T,
         network_color = 'red', spatial_network_name = 'Delaunay_network',
         point_size = 3, cell_color = 'Celltype_updated', cell_color_code = c("grey","#990099","orange","#0000FF","black"),
         title = "Delanuay triangulation network")

ggsave(plot_neighbors_path, plt, device = "png", width = 50, height = 25, units = "cm")

# network-averaging: smoothens the gene expression matrix by averaging the expression within one cell by using the neighbours within the predefined spatial network. 
# Instead of using k neighbors , we used Delanuay triangulation with specific radius. So the smoothin is done for all cells within radius
# Here we are computing correlation of gene expression along all combinations of genes within their network
spatialCorGenes = detectSpatialCorFeats(giotto_obj,
                                        expression_values = "normalized",
                                        method = "network",
                                        cor_method = "spearman")

# Statistical framework to identify if pairs of genes (such as ligand-receptor combinations) are expressed at higher 
# levels than expected based on a reshuffled null distribution of gene expression values in cells 
# that are spatially in proximity to each other..
set.seed(1)
CCI = spatCellCellcom(gobject = giotto_obj,
                spatial_network_name = "Delaunay_network",
                cluster_column = "Celltype",
                random_iter = 100, # number of permutation
                feat_set_1 = LR_database$ligand,
                feat_set_2 = LR_database$receptor,
                min_observations = 0,
                set_seed = 1,
                adjust_target = "feats")

# select cell type combination of interest
giotto_out_CT1_CT2 = CCI %>%  
  filter(lig_cell_type == "CT1" & rec_cell_type == "CT2") %>% 
  mutate(ligand_receptor = gsub("-","_",LR_comb))

# Average the Giotto results for cases like L-R1-R2 (aggregating L-R1 and L-R2)
giotto_LRR_averaged_out = list()
for(entry in unique(LR_database$ligand_receptor))
{
  name = str_split(entry,"_") %>% 
    unlist
  
  # generate all possible combinations between L and R (for matching with giotto output)
  n =  name %>% 
      expand.grid(.,.) %>% 
      filter(Var1!=Var2) %>% 
      apply(.,1,function(x) {str_c(x,collapse = "_")}) 
  
  # average spatial correlation between pairs of LR -> NOT USED FOR NOW
  spatial_correlation = map(giotto_out_CT1_CT2 %>% 
                              filter(ligand_receptor %in% n) %>% select(LR_comb) %>% unlist, function(x) {
    x  %<>% str_split(.,"-") %>% 
      unlist 
    
    x = spatialCorGenes$cor_DT %>% filter(feat_ID == x[1] & variable == x[2])
    
  }) %>% do.call(rbind.data.frame,.) %>% summarise(
    across(where(is.character), ~ first(.)), 
    across(where(is.factor), ~ first(.)), 
    across(where(is.numeric), ~ mean(., na.rm = TRUE))) 
  
  
  # Average over the results
  averaged_results = giotto_out_CT1_CT2 %>% 
    filter(ligand_receptor %in% n) %>%
    summarise(
      across(where(is.character), ~ first(.)),
      across(where(is.factor), ~ first(.)), 
      across(where(is.numeric), ~ mean(., na.rm = TRUE))) %>%
    mutate(spat_spearman_cor = spatial_correlation$spat_cor)
  
  giotto_LRR_averaged_out[[entry]] = averaged_results
} 
giotto_LRR_averaged_out = do.call(rbind.data.frame,giotto_LRR_averaged_out) %>%
  mutate(ligand_receptor = rownames(.))


# PI = log2fc * -log10(p.adj)
# log2fc > 0 means that A-B celltypes colocalize more often than by chance 
giotto_LRR_averaged_out %<>% select(ligand_receptor,p.adj,log2fc,PI,spat_spearman_cor) %>%
  mutate(significant = p.adj < 0.05,
         statistics = PI) %>%
  arrange(desc(PI))


  
# save data
write.table(giotto_LRR_averaged_out ,significant_interactions_path)