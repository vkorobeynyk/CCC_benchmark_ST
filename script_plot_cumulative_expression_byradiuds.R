
# Load package
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(sf)
  library(stringr)
  library(jsonlite)
  library(ggpubr)
  library(magrittr)
  library(purrr)
  source("scripts/helper_functions.R")
})

inflated_counts = read.csv("output/Visium_HD_HPC_semiSimulation_NB//inflated_normalized_counts_FC_1_FC_nSenderCells_1_FC_nReceiverCells_1_indexLR_1.tsv",sep="\t") %>% as.matrix
cellmetadata = read_json(path = "output/Visium_HD_HPC_semiSimulation_NB//simulated_cellmetadata_FC_1_FC_nSenderCells_1_FC_nReceiverCells_1_indexLR_1.json")

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

##########################################################################################
### Plot cumulative expression across distance for simulated LR and 1 significant pair ###
##########################################################################################

# select randomly same number of cells for every celltype
metadata = cellmetadata$metadata %>%
  group_by(Celltype) %>%
  sample_n(size = min(table(cellmetadata$metadata$Celltype))) %>%
  ungroup()
inflated_counts_subset = inflated_counts[,metadata$Cell_ID]

LR_database = read.table("data/LR_database.tsv")
genes = LR_database$ligand_receptor %>% str_split("_") %>% unlist %>% unique
genes = genes[genes %in% rownames(inflated_counts)]
plots_lst = list()
lst = list()
lst_anova = list()
lst_pvals = list()
for(gene in genes[1:2])
{
  # split the data into groups of 50 cells for each celltype to compute statistics in the end
  seeds = seq(1:20)
  for(seed in seeds)
  {
    x = cumulative_expression(gene = gene,seed = seed, inflated_counts = inflated_counts, cellmetadata = cellmetadata)
    plots_lst[[gene]][[paste0("seed_",seed)]] = x$plot
    lst[[gene]][[paste0("seed_",seed)]] = x$distance_curve_ct %>% mutate(seed = paste0("seed_",seed))
  }
  # compute statistics for each radius
  df = do.call(rbind.data.frame, lst[[gene]])
  
  library(broom)
  
  # split data by radius
  dfs = split(df, df$radius)
  
  # run ANOVA for each radius
  lst_pvals[[gene]] = map(dfs, ~ {
    x = t.test(.x %>% filter(sender_receiver == "CT1 CT2") %>% pull(total_expr), 
               .x %>% filter(sender_receiver != "CT1 CT2") %>% pull(total_expr), alternative = "greater")$p.value
  }) %>% as.data.frame %>% t
}

df = do.call(cbind.data.frame, lst_pvals)
df[is.na(df)] = 1
N_genes_CT1CT2_significant = rowSums(df < 0.05)
names(N_genes_CT1CT2_significant) = gsub("X", "", names(N_genes_CT1CT2_significant)) %>% as.integer
N_genes_CT1CT2_significant = data.frame(radius = names(N_genes_CT1CT2_significant),
                                        n_significant = N_genes_CT1CT2_significant)

barplot(N_genes_CT1CT2_significant$n_significant,
        names.arg = N_genes_CT1CT2_significant$radius,
        xlab = "radius",
        ylab = "count of significant results",
        main = "Significant t.test for CT1 CT2",
        col = "lightblue")

ggarrange(plotlist = list(plots_lst[[1]]$seed_1, 
                          plots_lst[[2]]$seed_1 , 
                          plots_lst[[3]]$seed_1,
                          plots_lst[[4]]$seed_1))
