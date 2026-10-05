# ============================================================================
# For every dataset / method / parameter combination (within one semi-
# simulation strategy - spatialScattering or spatialColocalization, see
# config.yaml), loads the method's significant_interactions output, checks
# whether the simulated Sender -> Receiver LR pair was recovered (and at
# what rank), and computes precision/recall/F1 treating that single
# simulated pair as ground truth. Results across all combinations are
# stored in a nested list and saved as output/{strategy}/final_scores.RDS
# ============================================================================

suppressMessages({
  library(dplyr)
  library(stringr)
  library(magrittr)
  library(glue)
})

# Fail early with a clear error if a required input/output is missing
if (is.null(snakemake@input[["significant_interactions"]]) | is.null(snakemake@output[["final_scores"]])) {
  stop("Argument_name needs to be specified, but is missing.n", call. = FALSE)
}

# ==============================================================================
# STEP 0: Snakemake I/O and params
# ==============================================================================

# output files
path_final_scores = snakemake@output[["final_scores"]]

# params
strategy = snakemake@wildcards[["strategy"]] # spatialScattering or spatialColocalization - see config.yaml

config = yaml::read_yaml("config.yaml")

methods = config$methods %>% unlist
datasets = config$datasets %>% unlist
indexLR_toSample = config$indexLR_toSample %>% unlist %>% as.double
radius_param_index = config$radius_param_index %>% unlist %>% as.double
PCE_Sender = config$semiSimulation$PCE_Sender %>% unlist %>% as.double
PCE_Receiver = config$semiSimulation$PCE_Receiver %>% unlist %>% as.double

# ==============================================================================
# STEP 1: loop over every dataset / indexLR / method / parameter combination
# and score each one against the known simulated ground truth
# ==============================================================================

statistics_lst = list()

for (dataset in datasets)
{
  for (indexLR in indexLR_toSample)
  {
    # load the simulated interaction for this indexLR (same indexLR always uses the same LR pair)
    simulated_interaction_file = list.files(paste0("output/", strategy, "/", dataset, "_semiSimulation_NB"), pattern = paste0("indexLR_", indexLR, ".RDS"))[1]
    simulated_interactions = readRDS(paste0("output/", strategy, "/", dataset, "_semiSimulation_NB/", simulated_interaction_file)) %>%
      .$Sender_Receiver
    
    for (method in methods)
    {
      # load all output files for every parameter combination, for this method
      data_dir = file.path("output", strategy, dataset, method)
      
      # build the file path for every (PCE_Sender, PCE_Receiver, radius_param_index) combination
      params_grid = expand.grid(indexLR = indexLR, PCE_Sender = PCE_Sender, PCE_Receiver = PCE_Receiver, radius_param_index = radius_param_index) %>%
        mutate(file_path = glue("{data_dir}/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR}.tsv"))
      
      for (file in params_grid$file_path)
      {
        x = params_grid %>% filter(file_path == file)
        naming = paste0("PCE_Sender_", x["PCE_Sender"], "_PCE_Receiver_", x["PCE_Receiver"], "_l_", x["radius_param_index"], "_indexLR_", x["indexLR"])
        
        all_interactions = read.table(file, header = T, sep = "\t")
        
        if (method == "mistyR") {
          all_interactions$importances = abs(all_interactions$importances) # we include both positive and negative correlations
        }
        
        significant_interactions = all_interactions %>% filter(as.logical(significant))
        
        # ------------------------------------------------------------------
        # rank of the simulated LR pair among this method's output
        # ------------------------------------------------------------------
        l_r = simulated_interactions %>% str_split("_") %>% unlist %>% str_c(., collapse = "_")
        
        l_r_all = grepl(l_r, all_interactions$ligand_receptor)
        if (nrow(significant_interactions) != 0)
        {
          l_r_significant = any(grepl(l_r, significant_interactions$ligand_receptor))
        } else {
          l_r_significant = FALSE
        }
        
        if (any(l_r_all)) {
          n = which(l_r_all) # some methods detect the L-R pair without detecting a specific subunit; in that case, take the best rank
          
          if (method == "cellchat" | method == "mistyR") {
            # methods where a higher statistic is more relevant -> rank in descending order
            rank = rank(all_interactions$statistics, ties.method = c("average"))[nrow(all_interactions) - n + 1] %>% min
          } else {
            rank = rank(all_interactions$statistics, ties.method = c("average"))[n] %>% min
          }
          
        } else {
          rank = NA # the inflated LR pair wasn't found in the output at all
        }
        
        # ------------------------------------------------------------------
        # precision / recall / F1, treating the single simulated LR pair as
        # the only ground-truth positive
        # ------------------------------------------------------------------
        TP = if (l_r_significant) {1} else {0}
        
        FP = if (l_r_significant) {
          nrow(significant_interactions) - 1 # inflated LR was significant -> FP = all other significant LR pairs found
        } else {
          nrow(significant_interactions) # inflated LR wasn't significant -> FP = every significant LR pair found
        }
        
        FN = if (l_r_significant) {0} else {1} # inflated LR was significant -> no false negative
        
        if (TP == 0) {
          precision = 0; recall = 0
        } else {
          precision = TP / (TP + FP); recall = TP / (TP + FN)
        }
        
        f1score = 2 * precision * recall / (precision + recall) %>% round(3)
        if (is.nan(f1score)) { f1score = 0 }
        
        # store the results for this combination
        df_statistics = data.frame(precision = precision %>% round(3), recall = recall %>% round(3), f1score = f1score %>% round(3),
                                   TP = TP, FN = FN, FP = FP, rank = rank, N = nrow(all_interactions),
                                   method = method, indexLR = indexLR,
                                   PCE_Sender = x["PCE_Sender"],
                                   PCE_Receiver = x["PCE_Receiver"],
                                   radius_param_index = x["radius_param_index"],
                                   ratio_Receiver_seen_byMethod = unique(all_interactions$ratio_Receiver_seen_byMethod),
                                   average_cells_perSender_seen_byMethod = unique(all_interactions$average_cells_perSender_seen_byMethod))
        
        print(paste("precision:", df_statistics$precision, "recall:", df_statistics$recall, "f1score:", df_statistics$f1score))
        
        statistics_lst[[dataset]][[method]][[naming]] = df_statistics
      }
    }
  }
}

# ==============================================================================
# STEP 2: save output
# ==============================================================================

saveRDS(statistics_lst, path_final_scores)
