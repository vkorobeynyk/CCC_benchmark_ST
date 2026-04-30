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
  library(ggpubr)
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
indexLR = snakemake@params["indexLR"] %>% as.numeric
FC_nReceiverCells = snakemake@params["FC_nReceiverCells"] %>% as.numeric
FC_nSenderCells = snakemake@params["FC_nSenderCells"] %>% as.numeric
radius = config[["l_param"]][["giotto"]][[dataset]][l_index] %>% unlist

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)
#############
# load data #
#############

inflated_counts = read.csv(normalized_counts_path,sep="\t") %>% as.matrix
cellmetadata = read_json(path = cellmetadata_path)

'
inflated_counts = read.csv("output/Visium_HD_HPC_semiSimulation_NB//inflated_normalized_counts_FC_1_FC_nSenderCells_4_FC_nReceiverCells_4_indexLR_1.tsv",sep="\t") %>% as.matrix
cellmetadata = read_json(path = "output/Visium_HD_HPC_semiSimulation_NB//simulated_cellmetadata_FC_1_FC_nSenderCells_4_FC_nReceiverCells_4_indexLR_1.json")
LR_database = read.table("data/LR_database.tsv", row.names = 1)
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

# CURRENTLY NOT USED
# network-averaging: smoothens the gene expression matrix by averaging the expression within one cell by using the neighbours within the predefined spatial network. 
# Instead of using k neighbors , we used Delanuay triangulation with specific radius. So the smoothin is done for all cells within radius
# Here we are computing correlation of gene expression along all combinations of genes within their network
spatialCorGenes = detectSpatialCorFeats(giotto_obj,
                                        expression_values = "normalized",
                                        spatial_network_name = "Delaunay_network",
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

##########################################################################################
### Plot cumulative expression across distance for simulated LR and 1 significant pair ###
##########################################################################################
if(FC_nReceiverCells == 1 & FC_nSenderCells == 1)
{
  simulated_genes = readRDS(paste0("output/", dataset, "_semiSimulation_NB/simulated_interactions_FC_1_FC_nSenderCells_1_FC_nReceiverCells_1_indexLR_" ,indexLR , ".RDS")) %>% unlist
  simulated_genes = c(simulated_genes , giotto_LRR_averaged_out %>% arrange(p.adj) %>% slice(1) %>% pull(ligand_receptor)) %>% str_split("_") %>% unlist
  plots_lst = list()
  for(gene in simulated_genes)
  {
    net = getSpatialNetwork(giotto_obj)
    df = net@networkDT
    
    # named vector for direct lookup
    celltype_vec = setNames(cellmetadata$metadata$Celltype, cellmetadata$metadata$cell_ID)
    
    df$ctype_from = celltype_vec[df$from]
    df$ctype_to   = celltype_vec[df$to]
    
    gene_expr = inflated_counts[gene, ]
    
    df$expr_from = gene_expr[df$from]
    df$expr_to   = gene_expr[df$to]
    df$expr_sum = df$expr_from + df$expr_to
    
    celltypes = unique(cellmetadata$metadata$Celltype)
    
    # you can change bin_size if needed
    bin_size = 10
    max_dist = max(df$distance)
    
    dist_bins = seq(0, max_dist, by = bin_size)
    
    library(dplyr)
    
    curve_list = lapply(celltypes, function(ct) {
      
      # Select edges where either side has the cell type of interest
      edges_ct = df %>% 
        filter(ctype_from == ct | ctype_to == ct)
      
      # Compute cumulative expression for each radius
      data.frame(
        celltype = ct,
        radius = dist_bins,
        total_expr = sapply(dist_bins, function(d)
          sum(edges_ct$expr_sum[edges_ct$distance <= d], na.rm = TRUE)
        )
      )
    })
    
    distance_curve_ct = bind_rows(curve_list)
    
    library(ggplot2)
    
    p1 = ggplot(distance_curve_ct, aes(x = radius, y = total_expr, color = celltype)) +
      geom_line(size = 1.2) +
      geom_point() +
      theme_classic(base_size = 14) +
      labs(
        x = "Distance radius",
        y = paste("Summed expression of", gene),
        title = paste("Cumulative gene expression by radius by cell type:", gene)
      ) +
      scale_color_brewer(palette = "Dark2") 
    plots_lst[[gene]] = p1
  }
  ggsave(paste0("output/Visium_HD_HPC/giotto/plot_cumulative_expressionbyRadius_l_", l_index), ggarrange(plotlist = plots_lst), device = "png", width = 50, height = 25, units = "cm")
}
  
# save data
write.table(giotto_LRR_averaged_out ,significant_interactions_path)