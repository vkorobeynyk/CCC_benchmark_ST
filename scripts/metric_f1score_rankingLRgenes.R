suppressMessages({
  library(dplyr)
  library(stringr)
  library(magrittr)
})

# An useful error if the argument is missing
if (is.null(snakemake@input[["significant_interactions"]]) | is.null(snakemake@input[["simulated_interactions"]]) | is.null(snakemake@output[["CT_statistics"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# OUTPUT FILES
path_CT_statistics <- snakemake@output[["CT_statistics"]]
path_ranking_LRgenes <- snakemake@output[["ranking_LRgenes"]]

# INPUT FILES
path_significant_interactions <- snakemake@input[["significant_interactions"]]
path_simulated_interactions <- snakemake@input[["simulated_interactions"]]

# PARAMS
method = snakemake@params[["method"]]
########################################
##### Loading and processing files #####
########################################

# load significant interactions
all_interactions = read.table(path_significant_interactions, header = T)
all_interactions$ligand_receptor = gsub("-","_",all_interactions$ligand_receptor)
# giotto allows to compute A-B neighbors that are spatially enriched or depleted. To be fair with other methods, we are only interested in enriched
if(method == "giotto") { 
  all_interactions %<>% filter(log2fc > 0)  # select only significant/positively enriched interactions compared to null
} else if(method == "mistyR") { 
  all_interactions$importances = abs(all_interactions$importances) # we are including positive and negative correlations
}
# select significant interactions}
significant_interactions = all_interactions %>% filter(as.logical(significant))

# load simulated interactions (TRUE)
simulated_interactions = readRDS(path_simulated_interactions) %>% 
  .$CT1_CT2

#############################################
##### Generating a ranking for LR genes #####
#############################################

lst_score_perCTCT = list()

# which index is the LR inflated

l_r = simulated_interactions %>% str_split("_") %>% unlist %>% str_c(., collapse = "_")

# only compute n when there are significant interactions
l_r_all = grepl(l_r, all_interactions$ligand_receptor)
if(nrow(significant_interactions) != 0) 
{
  l_r_significant = any(grepl(l_r, significant_interactions$ligand_receptor))
} else {l_r_significant = FALSE}

# This is needed as some method may have same statistic for different LR genes
if(any(l_r_all)) { 
  n = which(l_r_all) # there are cases where them method detect L-R without detecting the subunit. In this case, take the better rank
  if(method == "giotto" | method == "cellchat" | method == "mistyR") # methods that have statistics where the highest is the most relevant
  {
    rank = rank(all_interactions$statistics,
                ties.method = c("average"))[nrow(all_interactions)-n+1] %>% min
  } else {
    rank = rank(all_interactions$statistics,
                ties.method = c("average"))[n] %>% min
  }
  
} else{rank = NA} # if our inflated LR pair was not found in the output, rank NA


####################################################
##### Generating precision,recall and f1scores #####
####################################################

lst_score_perCTCT = list()

# compute statistics -> This is suitable for cases only we are inflating 1 ligand receptor pair
TP = if(l_r_significant) {1} else {0}

FP = if(l_r_significant) { 
  nrow(significant_interactions)-1 # if inflated LR was found in significant, FP = all significant LR found - 1
  
} else {nrow(significant_interactions)} # if inflated LR was not found in significant, FP = all significant LR found

FN = if(l_r_significant) {0} else {1} # if inflated LR was found in significatn, FN = 0 

if(TP == 0) {precision = 0 ; recall = 0} else {precision = TP / (TP + FP) ; recall = TP / (TP + FN)}
f1score = 2 * precision * recall / (precision + recall) %>% round(3)
if(is.nan(f1score)) {f1score = 0}

# Create df to store average results
df_statistics = data.frame(precision = precision %>% round(3), recall = recall %>% round(3), f1score = f1score %>% round(3), TP = TP, FN = FN, FP = FP)

# print all scores
print(paste("precision:", df_statistics$precision, "recall:" ,df_statistics$recall, "f1score:", df_statistics$f1score ))

# in case f1score is NaN
if(is.nan(df_statistics$f1score)) {df_statistics$f1score = 0}

######################
##### Save files #####
######################

write.table(df_statistics, path_CT_statistics, row.names = F)
write.table(rank, path_ranking_LRgenes, col.names = F, row.names = F)