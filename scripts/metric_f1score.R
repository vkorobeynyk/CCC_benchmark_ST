# Load package
library(dplyr)
library(jsonlite)
library(stringr)

# An useful error if the argument is missing
if (is.null(snakemake@input[["significant_interactions"]]) | is.null(snakemake@input[["simulated_interactions"]]) | is.null(snakemake@output[["CT_statistics"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# OUTPUT FILES
path_CT_statistics <- snakemake@output[["CT_statistics"]]

# INPUT FILES
path_significant_interactions <- snakemake@input[["significant_interactions"]]
path_simulated_interactions <- snakemake@input[["simulated_interactions"]]

########################################
##### Loading and processing files #####
########################################

# load significant interactions
significant_interactions = read.table(path_significant_interactions)

# load simulated interactions (TRUE)
simulated_interactions = readRDS(path_simulated_interactions) %>% .$CT1_CT2

####################################################
##### Generating precision,recall and f1scores #####
####################################################

lst_score_perCTCT = list()

simulated_interactions = simulated_interactions[!grepl("subunit", simulated_interactions)]

# compute statistics
TP = intersect(simulated_interactions, significant_interactions$ligand_receptor) %>% length
#TN = 
FP = setdiff(significant_interactions$ligand_receptor , simulated_interactions) %>% length
FN = setdiff(simulated_interactions , significant_interactions$ligand_receptor) %>% length

if(TP == 0) {precision = 0 ; recall = 0} else {precision = TP / (TP + FP) ; recall = TP / (TP + FN)}
f1score = 2 * precision * recall / (precision + recall) %>% round(2)

# Create df to store average results
df_statistics = data.frame(precision = precision %>% round(2), recall = recall %>% round(2), f1score = f1score %>% round(2))

# print all scores
print(paste("precision:", df_statistics$precision, "recall:" ,df_statistics$recall, "f1score:", df_statistics$f1score ))

# in case f1score is NaN
if(is.nan(df_statistics$f1score)) {df_statistics$f1score = 0}

######################
##### Save files #####
######################

write.csv(df_statistics, path_CT_statistics, row.names = F)