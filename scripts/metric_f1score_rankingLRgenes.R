# Load package
library(dplyr)
library(stringr)
library(magrittr)

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

########################################
##### Loading and processing files #####
########################################

# load significant interactions
all_interactions = read.table(path_significant_interactions, header = T)
all_interactions$ligand_receptor = gsub("-","_",all_interactions$ligand_receptor)
significant_interactions = all_interactions[all_interactions$significant %>% as.logical, ] # select only significant interactions

# load simulated interactions (TRUE)
simulated_interactions = readRDS(path_simulated_interactions) %>% 
  .$CT1_CT2

#############################################
##### Generating a ranking for LR genes #####
#############################################

lst_score_perCTCT = list()

# which index is the LR inflated

l_r = simulated_interactions %>% str_split("_") %>% unlist

n = grep(str_c(l_r, collapse = "_"), all_interactions$ligand_receptor)

if(length(n) != 0) { # if our inflated LR pair was not found in the output, rank NA
  n = min(n) # there are cases where them method detect L-R without detecting the subunit. In this case, take the better rank
  rank = rank(all_interactions$statistics,
              ties.method = c("average"))[n]
} else{rank = NA}


####################################################
##### Generating precision,recall and f1scores #####
####################################################

lst_score_perCTCT = list()

# compute statistics -> This is suitable for cases only we are inflating 1 ligand receptor pair
TP = if(length(n) != 0) {1} else {0}
FP = if(length(n) != 0) {
  if(nrow(significant_interactions) != 0) {nrow(significant_interactions)-1} else {nrow(significant_interactions)} # case when there are no significant interactions
} else {0}
FN = if(length(n) != 0) {0} else {1}

if(TP == 0) {precision = 0 ; recall = 0} else {precision = TP / (TP + FP) ; recall = TP / (TP + FN)}
f1score = 2 * precision * recall / (precision + recall) %>% round(3)

# Create df to store average results
df_statistics = data.frame(precision = precision %>% round(3), recall = recall %>% round(3), f1score = f1score %>% round(3), TP = TP, FN = FN, FP = FP)

# print all scores
print(paste("precision:", df_statistics$precision, "recall:" ,df_statistics$recall, "f1score:", df_statistics$f1score ))

# in case f1score is NaN
if(is.nan(df_statistics$f1score)) {df_statistics$f1score = 0}

######################
##### Save files #####
######################

write.csv(df_statistics, path_CT_statistics, row.names = F)
write.table(rank, path_ranking_LRgenes, col.names = F, row.names = F)