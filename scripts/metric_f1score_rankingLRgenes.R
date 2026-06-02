suppressMessages({
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(glue)
})

# An useful error if the argument is missing
if (is.null(snakemake@input[["significant_interactions"]]) |  is.null(snakemake@output[["final_scores"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# OUTPUT FILES
path_final_scores = snakemake@output[["final_scores"]]

##########
# PARAMS #
##########

config = yaml::read_yaml("config.yaml")

# Setting parameters
methods = config$methods %>% unlist
datasets = config$datasets %>% unlist
indexLR_toSample = config$indexLR_toSample %>% unlist %>% as.double
l_param_index = config$l_param_index %>% unlist %>% as.double
FC = config$semiSimulation$FC %>% unlist %>% as.double
FC_nSenderCells = config$semiSimulation$FC_nSenderCells %>% unlist %>% as.double
FC_nReceiverCells = config$semiSimulation$FC_nReceiverCells %>% unlist %>% as.double

########################################
##### Loading and processing files #####
########################################
statistics_lst = list()
for(dataset in datasets)
{
  for(indexLR in indexLR_toSample)
  {
    
    # load simulated interactions -> same indexLR has same interaction
    simulated_interaction_file = list.files(paste0("output/",dataset,"_semiSimulation_NB"), pattern = paste0("indexLR_",indexLR,".RDS"))[1]
    simulated_interactions = readRDS(paste0("output/",dataset,"_semiSimulation_NB/",simulated_interaction_file )) %>% 
      .$CT1_CT2
    
    for(method in methods)
    {
      # load all output files for every combination of parameter for specific method
      data_dir =  file.path("output/",dataset,method)
      
      # iterate over every combination of parameters and retrieve statistics
      params_grid = expand.grid(indexLR = indexLR, FC_nSenderCells = FC_nSenderCells, FC_nReceiverCells = FC_nReceiverCells,l_param_index = l_param_index) %>%
        mutate(file_path = glue("{data_dir}/significant_interactions_FC_1_FC_nSenderCells_{FC_nSenderCells}_FC_nReceiverCells_{FC_nReceiverCells}_l_{l_param_index}_indexLR_{indexLR}.tsv"))
      
      # iterate over the grid of parameters
      for(file in params_grid$file_path)
      {
        x = params_grid %>% filter(file_path == file)
        
        naming = paste0("FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"], "_l_", x["l_param_index"], "_indexLR_",x["indexLR"])
        
        # load correct count file depending on the params_grid
        all_interactions = read.table(file, header = T, sep ="\t") 
        
        # giotto allows to compute A-B neighbors that are spatially enriched or depleted. To be fair with other methods, we are only interested in enriched
        if(method == "giotto") { 
          all_interactions %<>% filter(log2fc > 0)  # select only significant/positively enriched interactions compared to null
        } else if(method == "mistyR") { 
          all_interactions$importances = abs(all_interactions$importances) # we are including positive and negative correlations
        }
        
        # select significant interactions}
        significant_interactions = all_interactions %>% filter(as.logical(significant))
        
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
        df_statistics = data.frame(precision = precision %>% round(3), recall = recall %>% round(3), f1score = f1score %>% round(3), 
                                   TP = TP, FN = FN, FP = FP, rank = rank, 
                                   method = method, indexLR = indexLR,
                                   FC_nSenderCells = x["FC_nSenderCells"], 
                                   FC_nReceiverCells = x["FC_nReceiverCells"], 
                                   l_param_index = x["l_param_index"],
                                   ratio_CT2_seen_byMethod = unique(all_interactions$ratio_CT2_seen_byMethod),
                                   average_cells_perCT1_seen_byMethod = unique(all_interactions$average_cells_perCT1_seen_byMethod))
        
        # print all scores
        print(paste("precision:", df_statistics$precision, "recall:" ,df_statistics$recall, "f1score:", df_statistics$f1score ))
        
        # in case f1score is NaN
        if(is.nan(df_statistics$f1score)) {df_statistics$f1score = 0}
        
        statistics_lst[[dataset]][[method]][[naming]] = df_statistics
      }
    }
  }
}

######################
##### Save files #####
######################

saveRDS(statistics_lst, path_final_scores)
