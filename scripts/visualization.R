suppressMessages({
  # Load package
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  library(ggrepel)
  library(reshape2)
  library(purrr)
  library(magrittr)
  library(plyr)
  library(scales)
  library(edgeR)
  library(jsonlite)
  library(SpatialExperiment)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(patchwork)
  library(ggspavis)
  library(sf)
  library(fmsb)
  library(ggridges)
  library(viridis)
  source("scripts/helper_functions.R")
})

# Call the argument
path_config.yaml = opt$config.yaml
path_output_dir = "output"
path_results_dir = "output/results"
dir.create(path_results_dir)

########################################
##### Loading and processing files #####
########################################

metric_results = readRDS("output/final_scores.RDS")

config = yaml::read_yaml("config.yaml")

# Setting parameters
methods = config$methods %>% unlist
datasets = config$datasets %>% unlist
indexLR = config$indexLR_toSample %>% unlist %>% as.double
l_param_index = config$l_param_index %>% unlist %>% as.double
FC = config$semiSimulation$FC %>% unlist %>% as.double
FC_nSenderCells = config$semiSimulation$FC_nSenderCells %>% unlist %>% as.double
FC_nReceiverCells = config$semiSimulation$FC_nReceiverCells %>% unlist %>% as.double

############################
##### Diagnostic plots #####
############################
master_lst_diagnosticPlots = list()
diagnostic_gene_densityPlots = list()
diagnostic_plots_realFC = list()
diagnostic_plots_perCT = list()
gene_metadata_lst = list()
for(dataset in datasets)
{
  # Select inflated count files
  file_path = file.path(path_output_dir,paste0(dataset,"_semiSimulation_NB"))
  inflated_counts_files = file_path %>% list.files(., pattern = "inflated_counts")
  
  # Select all simulated interactions files
  simulated_interactions_files = list.files(file_path, pattern = "simulated_interactions")
  
  # Select simulated cellmetadata files
  cellmetadata_path = file.path(path_output_dir,paste0(dataset,"_semiSimulation_NB"))
  simulated_cellmetadata_files = list.files(cellmetadata_path, pattern = "simulated_cellmetadata")
  
  #### Read original counts for visium datasets
  original_counts = read.table(file.path("data/processed",dataset, paste0("processed_counts_",dataset,".tsv")))
  
  # generate the grid of parameters used for naming the list to generate outputs
  # I am not include FC param as it will be deprecated but it doesnt matter if I add it here
  params_grid = expand.grid(indexLR = indexLR, FC_nSenderCells = FC_nSenderCells, FC_nReceiverCells = FC_nReceiverCells)
  
  # iterate over the grid of parameters
  for(i in 1:nrow(params_grid))
  {
    x = params_grid[i,] %>% as.numeric ;  names(x) = c("indexLR","FC_nSenderCells","FC_nReceiverCells")
    
    naming = paste0("FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"], "_indexLR_",x["indexLR"])
    
    # load correct count file depending on the params_grid
    inflated_counts_file = inflated_counts_files %>% magrittr::extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_indexLR_",x["indexLR"] , ".tsv$"), inflated_counts_files))
    inflated_counts = read.table(file.path(file_path,inflated_counts_file)) 
    
    # load correct simulated_cellmetadata file depending on the params_grid
    simulated_cellmetadata_file = simulated_cellmetadata_files %>% magrittr::extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_indexLR_",x["indexLR"] , ".json$"), simulated_cellmetadata_files))
    master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]] = read_json(file.path(file_path,simulated_cellmetadata_file)) %>% convert_json_to_df
    
    # load correct simulated interactions file depending on the params_grid
    simulated_interactions_file = simulated_interactions_files %>% magrittr::extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_indexLR_",x["indexLR"] , ".RDS$"), simulated_interactions_files))
    master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] = readRDS(file.path(file_path,simulated_interactions_file))
    
    master_lst_diagnosticPlots[[naming]][["CT1"]] = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]$CT1_CT2$ligand
    master_lst_diagnosticPlots[[naming]][["CT2"]] = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]$CT1_CT2$receptor
    
    #############################################################################
    ##### Generate density plots of simulated genes before/after simulation #####
    
    for(CT in c("CT1","CT2"))
    {
      # check if cellnames are ordered
      stopifnot(master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$cell_ID == colnames(inflated_counts))
      
      set.seed(1)
      gene_names = master_lst_diagnosticPlots[[naming]][[CT]] %>% unlist
      
      # in case there are subunits, separate the genes
      if(grepl("_", gene_names)) {
        gene_names = str_split(gene_names, "_") %>% unlist
      } 
      
      for(gene in gene_names)
      {
        #####################################################################################
        ##### cumulative expression for cells where we added signal and all other cells #####
        
        m = master_lst_diagnosticPlots[[naming]]$simulated_cellmetadata$metadata
        
        gene_metadata_lst[[dataset]][[naming]][[CT]][[gene]] = data.frame(
          expression = inflated_counts[gene, ] %>% as.vector %>% unlist,
          celltype = m$Celltype 
        ) %>%
          group_by(celltype) %>%
          dplyr::summarise(
            # A cell expresses the gene if count > 0
            num_expressing = sum(expression > 0),
            total_cells = n(),
            indexLR = x["indexLR"]
          ) %>% filter(celltype == CT)
        
        # -----------------------------------
        
        tmp_metadata = master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$metadata
        # create df with counts of inflated and original count matrices
        # since we have CT1 and CT1_signalAdded, use grepl to find both
        tmp_df = data.frame(inflated_counts = inflated_counts[gene,grepl(CT, tmp_metadata$Celltype)] %>% as.numeric,
                            original_counts = original_counts[gene,grepl(CT, tmp_metadata$Celltype)] %>% as.numeric) %>% melt
        
        # calculate mean value for both counts
        mu = ddply(tmp_df, "variable", summarise, grp.mean=mean(value))
        
        p = ggplot(tmp_df, aes(x=value, color=variable)) +
          geom_density()+
          geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),
                     linetype="dashed") +
          ggtitle(paste0("density of counts | lines - mean | FC_nSenderCells:",x["FC_nSenderCells"], "FC_nReceiverCells " ,x["FC_nReceiverCells"] , "  gene:", gene)) +
          theme(plot.title = element_text(size = 7, face = "bold"),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.title=element_blank(),
                axis.text=element_text(size=6)) 
        
        diagnostic_gene_densityPlots[[dataset]][[naming]][[CT]][[gene]] = p
      }
      rm(tmp_metadata)
    }
    
    ##########################
    ##### Save AvelogCPM #####
    # This is for "compute_diagnostic_plots" function in order to save memory instead of saving entire inflated counts
    master_lst_diagnosticPlots[[naming]][["inflated_counts_aveLogCPM"]] = inflated_counts[,master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$metadata$Celltype %in% c("CT1", "CT2")] %>% aveLogCPM
    names(master_lst_diagnosticPlots[[naming]][["inflated_counts_aveLogCPM"]]) = rownames(inflated_counts)
  }
  
  ######################################################
  ##### Plot avelogCPM before and after simulation #####
  
  # filter original avelogcpm counts to contain same genes as the inflated count matrices
  original_counts = original_counts[which(rownames(original_counts) %in% rownames(inflated_counts)),]
  
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, indexLR = indexLR,
                                                               FC_nSenderCells = FC_nSenderCells, FC_nReceiverCells = FC_nReceiverCells, dataset = dataset, CT_toPlot = c("CT1","CT2")) # CT_toPlot has to be same as above
}

saveRDS(master_lst_diagnosticPlots, "output/results/master_lst_diagnosticPlots.RDS")
saveRDS(gene_metadata_lst , "output/results/gene_metadata_lst.RDS")

######################
##### Save plots #####
for(dataset in datasets)
{
  pdf(file.path(path_results_dir ,paste0(dataset, "_diagnostic_plots.pdf")), width = 12, height = 7)
  do.call(ggarrange,c(diagnostic_plots_perCT[[dataset]]$avelogcpm, common.legend = TRUE)) %>% print
  
  #x = 1:length(diagnostic_plots_realFC[[dataset]])
  #sapply(x, function(x) {do.call(ggarrange,diagnostic_plots_realFC[[dataset]][[x]]) %>% print}) %>% print
  #diagnostic_plots_MeanVar[[dataset]] %>% print
  
  # generate density plots
  N = length(diagnostic_gene_densityPlots[[dataset]])
  indices = round(seq(1, N, length.out = 4))
  
  do.call(ggarrange,c(diagnostic_gene_densityPlots[[dataset]][indices][[1]][["CT1"]],  
                      diagnostic_gene_densityPlots[[dataset]][indices][[2]][["CT2"]],
                      diagnostic_gene_densityPlots[[dataset]][indices][[3]][["CT1"]],  
                      diagnostic_gene_densityPlots[[dataset]][indices][[4]][["CT2"]], common.legend = TRUE)) %>% print
  
  #barplot(N_genes_CT1CT2_significant$n_significant,
  #        names.arg = N_genes_CT1CT2_significant$radius,
  #        xlab = "radius",
  #        ylab = "count of significant results",
  #        main = "Amount of Receptors for CT1-CT2 whose cumulative expression is statistically significant for each radius",
  #        col = "lightblue") %>% print
  
  #ggarrange(plotlist = list(cumulative_expression_acrossRadius[[1]]$seed_1, 
  #                          cumulative_expression_acrossRadius[[2]]$seed_5 , 
  #                          cumulative_expression_acrossRadius[[3]]$seed_8,
  #                          cumulative_expression_acrossRadius[[4]]$seed_15)) %>% print
  dev.off()
}

##################################
##### precision/recall plots #####
##################################

# plot TPR/sensitivity/recall
statistics_results_lst_recallprecision = list()
master_lst_precision_recall = list()
simulated_interactions_lst = list()
for(dataset in datasets)
{
  for(method in methods)
  {
    # Load files for upset plot
    significant_interactions_files = file.path(path_output_dir,dataset,method) %>% 
      list.files(., pattern = "significant_interactions")
    # Load files for upset plot
    simulated_interactions_files = file.path(path_output_dir,paste0(dataset,"_semiSimulation_NB")) %>% 
      list.files(., pattern = "simulated_interactions")
    
    for(file in significant_interactions_files)
    {
      # retrieve all significant interactions per method
      master_lst_precision_recall[["significant_interactions"]][[method]][[file]] = read.table(file.path(path_output_dir,dataset,method,file),header = T, sep = "\t") %>% 
        mutate(significant = as.logical(significant)) %>%
        filter(significant) %>% 
        select(ligand_receptor) %>% 
        unlist %>%
        unname
      
      # retrieve all significant among the simulated interactions per method
      indexLR = str_extract(file, "(?<=indexLR_)\\d+")
      
      n = simulated_interactions_files %>% grep(paste0("indexLR_",indexLR),.) %>%
        magrittr::extract(1) # the simulated interaction only changes with indexLR, so take the first
      
      # this stores same things multiple times
      # can be improved but for now it doesnt take much space nor time
      simulated_interactions_lst[[method]][[file]] = readRDS(file.path(path_output_dir,paste0(dataset,"_semiSimulation_NB"),simulated_interactions_files[n])) %>%
        unlist %>%
        str_flatten("_")
      
      master_lst_precision_recall[["intersect_significant_simulated"]][[method]][[file]] = ifelse(simulated_interactions_lst[[method]][[file]] %in% master_lst_precision_recall[["significant_interactions"]][[method]][[file]],
                                                                                                  simulated_interactions_lst[[method]][[file]], 0)
    }
    
    
    ################################################################################################
    ##### Average the metrics for different indexLR but across same combination of parameters  #####
    
    statistics_results_lst_recallprecision_non_indexLRaveraged = metric_results %>% lapply(., function(dataset) {
      # iterate over each method and average metrics of indexLR parameter
      lapply(dataset, function(method) {
        tmp_df = method %>% do.call(rbind, .) %>%
          # Create the new column by removing the suffix
          # The regex "_indexLR_\\d+$" targets "_indexLR_" and all digits at the end
          mutate(filename = rownames(.))
      })
    })
    
    # extract rows that have the same parameters besides indexLR and then average the columns. We are averaging results of indexLR
    statistics_results_lst_recallprecision = metric_results %>% lapply(., function(dataset) {
      # iterate over each method and average metrics of indexLR parameter
      lapply(dataset, function(method) {
        tmp_df = method %>% do.call(rbind, .) %>%
          # Create the new column by removing the suffix
          # The regex "_indexLR_\\d+$" targets "_indexLR_" and all digits at the end
          mutate(file_without_indexLR = str_remove(rownames(.), "_indexLR_\\d+$")) %>%
          group_by(file_without_indexLR) %>%
          # average every numeric column
          dplyr::summarise(
            across(
              everything(), 
              ~ if(is.numeric(.x)) mean(.x, na.rm = TRUE) else .x[1]
            ), 
            .groups = "drop"
          )
      })
    })
  }
  
  data_across_methods = statistics_results_lst_recallprecision[[dataset]]
  ###########################################################################
  ##### Calculate % of cells from cellmetadata that we added signal to  #####
  
  # we dont care about l_param_index parameter (here we use it for matching with master list) nor indexLR as that doesnt affect amount of cells meaning all files in simulated_cellmetadata_file should contain same %
  # load correct simulated_cellmetadata file depending on the params_grid | since all datasets have same values for % of cells signal added, which simulated_cellmetadata to load is not importnat 
  
  params_grid = expand.grid(FC_nSenderCells = FC_nSenderCells, FC_nReceiverCells = FC_nReceiverCells, l_param_index = l_param_index)
  tmp_lst = list()
  for(row in 1:nrow(params_grid))
  {
    x = params_grid[row,] %>% as.numeric ;  names(x) = c("FC_nSenderCells","FC_nReceiverCells", "l_param_index")
    
    naming = paste0("FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_l_", x["l_param_index"])
    
    simulated_cellmetadata_file = simulated_cellmetadata_files %>% magrittr::extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"]), simulated_cellmetadata_files))
    
    y = read_json(file.path(path_output_dir,paste0(dataset, "_semiSimulation_NB"),simulated_cellmetadata_file[1])) %>% convert_json_to_df # take only 1 file -> to speed up computations
    y = y$metadata
    # calculate percentage of "CT1_signal_added" and "CT2_signal_added" cells
    pct_signal_added = y %>%
      filter(Celltype != "Other") %>% 
      group_by(Celltype) %>%
      dplyr::summarise(
        pct_signal_added = round(100 * mean(grepl("signalAdded", Celltype_updated))),
        amount_signal_added = sum(grepl("signalAdded", Celltype_updated))
      )
    
    tmp_lst[[naming]] = pct_signal_added
  }
  # Reorder the vector to match the order of df$match_id
  # This looks up the values in my_vector based on the ID names in the dataframe
  reordered_list = tmp_lst[data_across_methods[[1]]$file_without_indexLR]
  
  # make sure the entire reordered vector has same names as the dataframe column
  stopifnot(!any(reordered_list %>% names %>% is.na))
  
  #### Combine all methods together
  # for precision recall plots
  # add pct_signal_added variable
  
  data_across_methods = do.call(rbind, data_across_methods) %>%
    mutate(pct_CT1cells_expressing_L = 
             rep(reordered_list %>% do.call(rbind.data.frame,.) %>% filter(Celltype == "CT1") %>% pull(pct_signal_added), length(methods)),
           N_CT1cells_expressing_L = 
             rep(reordered_list %>% do.call(rbind.data.frame,.) %>% filter(Celltype == "CT1") %>% pull(amount_signal_added), length(methods)),
           pct_CT2cells_expressing_R = 
             rep(reordered_list %>% do.call(rbind.data.frame,.) %>% filter(Celltype == "CT2") %>% pull(pct_signal_added), length(methods)),
           N_CT2cells_expressing_R = 
             rep(reordered_list %>% do.call(rbind.data.frame,.) %>% filter(Celltype == "CT2") %>% pull(amount_signal_added), length(methods))
    )
  
  ##################################################
  ##### UpSet plots all interactions retrieved #####
  
  significant_interactions_retrieved = lapply(master_lst_precision_recall[["significant_interactions"]] , function(method) {
    # find all unique elements
    all_strings = unique(unlist(master_lst_precision_recall$significant_interactions))
    # For each element in the list, check if the global strings exist there
    binary_df = as.data.frame(lapply(method, function(x) {
      as.numeric(all_strings %in% x)
    }))
    
    data.frame(method = ifelse(rowSums(binary_df) > 0, 1,0), row.names = all_strings)
  }) %>% do.call(cbind,.)
  
  colnames(significant_interactions_retrieved) = names(master_lst_precision_recall$significant_interactions)
  
  # create upset plot
  m = make_comb_mat(significant_interactions_retrieved)
  upset1 = UpSet(
    m, 
    top_annotation = upset_top_annotation(
      m, 
      height = unit(12.5, "cm"), 
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 12),
      gp = gpar(fill = "steelblue")
    ),
    
    pt_size = unit(3, "mm"),        # Smaller dots
    lwd = 1,                        # Thinner lines 
    
    right_annotation = upset_right_annotation(
      m, 
      width = unit(3, "cm"),        # Keep set-size bars narrow
      gp = gpar(fill = "darkred")
    ),
    
    row_names_gp = gpar(fontsize = 12), # Smaller font
    row_gap = unit(0, "mm"),           # Removes extra space between rows
    
    column_title = "All significant interactions retrieved across all combination of parameters",
    column_title_gp = gpar(fontsize = 16, fontface = "bold")
  )
  
  ###################################################
  ##### UpSet plots simulated among significant #####
  intersect_significant_simulated_interactions = lapply(master_lst_precision_recall[["intersect_significant_simulated"]] , function(method) {
    # find all unique elements
    all_strings = unique(unlist(simulated_interactions_lst))
    # For each element in the list, check if the global strings exist there
    binary_df = as.data.frame(lapply(method, function(x) {
      as.numeric(all_strings %in% x)
    }))
    
    data.frame(method = ifelse(rowSums(binary_df) > 0, 1,0), row.names = all_strings)
  }) %>% do.call(cbind,.)
  
  # remove 0
  intersect_significant_simulated_interactions = intersect_significant_simulated_interactions[rownames(intersect_significant_simulated_interactions) != "0",]
  
  colnames(intersect_significant_simulated_interactions) = names(master_lst_precision_recall$intersect_significant_simulated)
  
  # create upset plot
  m = make_comb_mat(intersect_significant_simulated_interactions)
  upset2 = UpSet(
    m, 
    top_annotation = upset_top_annotation(
      m, 
      height = unit(12.5, "cm"), 
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 12),
      gp = gpar(fill = "steelblue")
    ),
    
    pt_size = unit(3, "mm"),        # Smaller dots
    lwd = 1,                        # Thinner lines 
    
    right_annotation = upset_right_annotation(
      m, 
      width = unit(3, "cm"),        # Keep set-size bars narrow
      gp = gpar(fill = "darkred")
    ),
    
    # --- TEXT & SPACING ---
    row_names_gp = gpar(fontsize = 12), # Smaller font
    row_gap = unit(0, "mm"),           # Removes extra space between rows
    
    column_title = "All significant interactions retrieved across all combination of parameters",
    column_title_gp = gpar(fontsize = 16, fontface = "bold")
  )
  
  ##########################################
  ##### Generate Precision recall plot #####
  
  tmp_df = data_across_methods %>%
    mutate(ratio_ReceiverSender = (FC_nReceiverCells/FC_nSenderCells) %>% log2)
  
  p1 = ggplot(tmp_df, aes(x = precision, y = recall)) +
    geom_density_2d(
      aes(fill = after_stat(level)), 
      contour_var = "ndensity",
      h = c(0.05,0.05) # fine tude KDE bandwidth as some methods have very low variance and density estimation doesnt work
    ) +
    facet_wrap(~ method, ncol = 4) +
    labs(title = "Independent 2D Density per method")
  
  ###################################
  ##### Generate f1 score plots #####
  
  tmp_df = data_across_methods %>% 
    mutate(AcountsSpatialDistance = ifelse(method %in% c("cellphonedb", "seurat_wilcoxon"), FALSE, TRUE)) # generate column showing if the method is spatial distance dependent or not
  
  d = tmp_df  %>%
    group_by(method,l_param_index ) %>%
    summarise_at(vars(f1score), mean) %>%
    mutate(f1score = round(f1score,3))
  
  heatmap1 = ggplot(d, aes(method, l_param_index, fill= f1score)) +
    geom_tile(color = "black") +
    geom_text(aes(label = f1score), color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("heatmap of f1score according to index l (the value of l param can be seen in config.yaml)")+
    scale_y_continuous(breaks=l_param_index) +
    ylab("Index of l parameter")
  
  p3 = ggplot(tmp_df, aes(x = pct_CT1cells_expressing_L, y = f1score, color = method, linetype = AcountsSpatialDistance)) + 
    geom_point(size = 0.25) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l_param_index)+
    ggtitle("F1score faceted by L")  +
    scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid"))+
    theme_bw() +
    theme(legend.key.width = unit(1.5, "cm")) + 
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
    scale_color_viridis_d(option = "D", name = "Methods") 
  
  p4 = ggplot(tmp_df, aes(x = N_CT1cells_expressing_L, y = f1score, color = method, linetype = AcountsSpatialDistance)) + 
    geom_point(size = 0.25) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l_param_index)+
    ggtitle("F1score faceted by L")  +
    scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid"))+
    theme_bw() +
    theme(legend.key.width = unit(1.5, "cm")) + 
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
    scale_color_viridis_d(option = "D", name = "Methods") 
  
  p5 = ggplot(tmp_df, aes(x = pct_CT2cells_expressing_R, y = f1score, color = method, linetype = AcountsSpatialDistance)) + 
    geom_point(size = 0.25) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l_param_index)+
    ggtitle("F1score faceted by L")  +
    scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid"))+
    theme_bw()
  
  p6 = ggplot(tmp_df, aes(x = N_CT2cells_expressing_R, y = f1score, color = method, linetype = AcountsSpatialDistance)) + 
    geom_point(size = 0.25) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l_param_index)+
    ggtitle("F1score faceted by L")  +
    scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
    theme_bw()
  
  ##########################################
  ##### Generate heatmap across params #####
  
  tmp_df = data_across_methods
  
  #### alternative faceted heatmap
  heatmap2_faceted_f1score = 
    ggplot(tmp_df, aes(x = factor(FC_nSenderCells), y = factor(FC_nReceiverCells), fill = f1score)) +
    geom_tile() +
    facet_grid(l_param_index ~ method, labeller = label_both) + 
    scale_fill_viridis_c(option = "magma",limits = c(0, 1)) + # 'magma' is great for seeing F1 hotspots
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey90"), 
      panel.spacing = unit(0.5, "lines") 
    ) +
    labs(
      title = "Comparing methods across radius and Fold-Change thresholds",
      x = "FC_nSenderCells",
      y = "FC_nReceiverCells",
      fill = "F1-Score"
    )
  
  heatmap2_faceted_recall = ggplot(tmp_df, aes(x = factor(FC_nSenderCells), y = factor(FC_nReceiverCells), fill = recall)) +
    geom_tile() +
    facet_grid(l_param_index ~ method, labeller = label_both) + 
    scale_fill_viridis_c(option = "magma",limits = c(0, 1)) + # 'magma' is great for seeing F1 hotspots
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey90"),
      panel.spacing = unit(0.5, "lines") 
    ) +
    labs(
      title = "Comparing methods across radius and Fold-Change thresholds",
      x = "FC_nSenderCells",
      y = "FC_nReceiverCells",
      fill = "Recall"
    )
  
  # Create a 'Combination' label for columns
  tmp_df %<>%
    mutate(Param_Comb = paste0("FC_S:", FC_nSenderCells, "\nFC_R:", FC_nReceiverCells, "\nl:",l_param_index)) %>%
    select(method, Param_Comb, recall) %>%
    pivot_wider(names_from = Param_Comb, values_from = recall) %>%
    tibble::column_to_rownames("method") %>%
    as.matrix()
  
  # Extract metadata for annotations (the labels at the top)
  # This keeps the FC and PCE values linked to the columns
  col_meta = data.frame(colnames(tmp_df)) %>%
    separate(1, into = c("FC_S", "FC_R","l"), sep = "\n") 
  
  col_meta$FC_S = factor(col_meta$FC_S, levels = c("FC_S:2","FC_S:4","FC_S:6","FC_S:8","FC_S:10"))
  col_meta$FC_R = factor(col_meta$FC_R, levels = c("FC_R:2","FC_R:4","FC_R:6","FC_R:8","FC_R:10"))
  col_meta$l = factor(col_meta$l, levels = unique(col_meta$l))
  
  # order
  #col_meta$FC_S = factor(col_meta$FC_S, levels = c("FC_S:0.3","FC_S:0.5","FC_S:0.8","FC_S:1","FC_S:1.2","FC_S:1.6","FC_S:2","FC_S:2.5"))
  #col_meta$nFC_R = factor(col_meta$nFC_R, levels =  c("FC_R:0.3","FC_R:0.5","FC_R:0.8","FC_R:1","FC_R:1.2","FC_R:1.6","FC_R:2","FC_R:2.5"))
  #col_meta$l = factor(col_meta$l, levels = unique(col_meta$l))
  
  # Get unique values
  unique_fcs = unique(col_meta$FC_S)
  unique_fcr = unique(col_meta$FC_R)
  unique_l = unique(col_meta$l)
  
  fcs_cols = setNames(brewer.pal(length(unique_fcs), "Set1"), unique_fcs)
  #fcs_cols = fcs_cols[1:2]
  fcr_cols = setNames(brewer.pal(length(unique_fcr), "Set2"), unique_fcr)
  #fcr_cols = fcr_cols[1:2]
  l_cols = setNames(brewer.pal(length(unique_l), "Set2"), unique_l)
  #l_cols = l_cols[1:2]
  
  top_ann = HeatmapAnnotation(
    FC_S = col_meta$FC_S,
    FC_R = col_meta$FC_R,
    col = list(
      FC_S = fcs_cols, 
      FC_R = fcr_cols, 
      l = l_cols
    )
  )
  
  
  # f1 score color gradient (0 to 1)
  col_f1 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  
  # Define your target "canvas" size for the heatmap body in mm
  target_width_mm = 400
  target_height_mm = 40
  
  # Calculate independent sizes
  n_cols = ncol(tmp_df)
  n_rows = nrow(tmp_df)
  
  dynamic_cell_width  = target_width_mm / n_cols
  dynamic_cell_height = target_height_mm / n_rows
  
  heatmap2 = Heatmap(tmp_df, 
                     name = "Recall",
                     top_annotation = top_ann,
                     column_split = col_meta$l,
                     cluster_rows = FALSE,
                     cluster_columns = FALSE, 
                     
                     # --- INDEPENDENT DYNAMIC SIZING ---
                     width = n_cols * unit(dynamic_cell_width, "mm"),
                     height = n_rows * unit(dynamic_cell_height, "mm"),
                     
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 8),
                     show_column_names = FALSE,
                     
                     # Adjust border thickness based on cell height to prevent "bleeding"
                     rect_gp = gpar(col = "white", lwd = if(dynamic_cell_height < 3) 0 else 1)
  ) %>% ComplexHeatmap::draw(., 
                             heatmap_legend_side = "right", 
                             annotation_legend_side = "bottom",
                             padding = unit(c(10, 2, 2, 20), "mm"), # Adds space at the bottom
                             merge_legends = FALSE)
  #############################################################
  ##### Plot how many overall cells are seen by method ########
  ##### This piece of code must be under the above heatmap ####
  ##### as it takes all the information from above ############
  
  tmp_df = data_across_methods
  
  tmp_df %<>%
    filter(!method %in% c("seurat_wilcoxon","cellphonedb")) %>% # remove spatial unaware methods
    mutate(average_cells_perCT1_seen_byMethod = as.numeric(average_cells_perCT1_seen_byMethod),
           Param_Comb = paste0("FC_S:", FC_nSenderCells, "\nFC_R:", FC_nReceiverCells, "\nl:",l_param_index)) %>%
    select(method, Param_Comb, average_cells_perCT1_seen_byMethod) %>%
    pivot_wider(names_from = Param_Comb, values_from = average_cells_perCT1_seen_byMethod) %>%
    tibble::column_to_rownames("method") %>%
    as.matrix()
  
  top_ann = HeatmapAnnotation(
    FC_S = col_meta$FC_S,
    FC_R = col_meta$FC_R,
    col = list(
      FC_S = fcs_cols, 
      FC_R = fcr_cols, 
      l = l_cols
    )
  )
  # Define your target "canvas" size for the heatmap body in mm
  target_width_mm = 300
  target_height_mm = 30
  
  # Calculate independent sizes
  n_cols = ncol(tmp_df)
  n_rows = nrow(tmp_df)
  
  dynamic_cell_width  = target_width_mm / n_cols
  dynamic_cell_height = target_height_mm / n_rows
  
  heatmap3 = Heatmap(tmp_df, 
                     name = "average_cells_perCT1_seen_byMethod",
                     top_annotation = top_ann,
                     column_split = col_meta$l,
                     width = n_cols * unit(dynamic_cell_width, "mm"),
                     height = n_rows * unit(dynamic_cell_height, "mm"),
                     cluster_rows = FALSE,
                     cluster_columns = FALSE, 
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 8),
                     show_column_names = FALSE,
                     rect_gp = gpar(col = "white", lwd = 1)) %>%
    ComplexHeatmap::draw(., 
                         heatmap_legend_side = "right", 
                         annotation_legend_side = "bottom",
                         # Use 'adjusted' which is the default valid argument
                         legend_grouping = "adjusted", 
                         ht_gap = unit(10, "mm"),
                         merge_legends = FALSE)
  
  #############################################################
  ##### Plot how many CT2 cells are captured by method ########
  ##### This piece of code must be under the above heatmap ####
  ##### as it takes all the information from above ############
  
  tmp_df = data_across_methods
  
  tmp_df %<>%
    mutate(Param_Comb = paste0("FC_S:", FC_nSenderCells, "\nFC_R:", FC_nReceiverCells, "\nl:",l_param_index)) %>%
    select(method, Param_Comb, ratio_CT2_seen_byMethod) %>%
    pivot_wider(names_from = Param_Comb, values_from = ratio_CT2_seen_byMethod) %>%
    filter(!method %in% c("seurat_wilcoxon","cellphonedb")) %>% # remove spatial unaware methods
    tibble::column_to_rownames("method") %>%
    as.matrix()
  
  heatmap4 = Heatmap(tmp_df, 
                     name = "ratio_CT2_seen_byMethod",
                     top_annotation = top_ann,
                     column_split = col_meta$l, 
                     cluster_rows = FALSE,
                     cluster_columns = FALSE, 
                     width = n_cols * unit(dynamic_cell_width, "mm"),
                     height = n_rows * unit(dynamic_cell_height, "mm"),
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 8),
                     show_column_names = FALSE,
                     rect_gp = gpar(col = "white", lwd = 1)) %>%
    ComplexHeatmap::draw(., 
                         heatmap_legend_side = "right", 
                         annotation_legend_side = "bottom",
                         # Use 'adjusted' which is the default valid argument
                         legend_grouping = "adjusted", 
                         ht_gap = unit(10, "mm"),
                         merge_legends = FALSE)
  
  ##########################################
  ##### Generate ranking LR genes plot #####
  tmp_df = data_across_methods %>% 
    filter(!method %in% c("seurat_wilcoxon", "cellphonedb")) %>%
    mutate(average_cells_perCT1_seen_byMethod = as.numeric(average_cells_perCT1_seen_byMethod) %>% as.integer)
  
  # Reshape data so FP and f1score are in a single column
  tmp_df_long <- tmp_df %>%
    pivot_longer(
      cols = c(FP, f1score), 
      names_to = "metric", 
      values_to = "value"
    ) %>%
    mutate(metric = case_when(
      metric == "FP" ~ "Number of Interactions Retrieved",
      metric == "f1score" ~ "F1-Score"
    ))
  
  N_interactions_f1score_plot = ggplot(tmp_df_long, aes(x = average_cells_perCT1_seen_byMethod, y = value, color = method)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_smooth(se = FALSE, span = 0.5, size = 0.8) + 
    facet_grid(metric ~ l_param_index, scales = "free_y") +
    scale_y_log10() +
    labs(
      x = "average_cells_perCT1_seen_byMethod",
      y = NULL, 
      title = "Performance Metrics Across Parameter Combinations"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  
  # generate histogram of number of NaN per method (NaN is when the method doesnt have the simulated ligand-receptor in its output)
  d = tmp_df %>%
    group_by(method) %>%
    dplyr::summarise(n_na = sum(is.na(rank)))
  
  p2 = ggplot(d, aes(x = method, y = n_na, fill = method)) +
    geom_col() +
    labs(
      x = "Method",
      y = "Count",
      title = "Amount of NaN - method didnt find the simulated LR pair"
    ) +
    scale_y_continuous(breaks = 0:max(d$n_na))
  
  # Ridge plot of rank according to l_param_index radius
  ridge_plot = ggplot(tmp_df, aes(x = rank, y = method, fill = method)) +
    geom_density_ridges() +
    scale_x_log10() + 
    theme_bw() +
    ggtitle("rank of how each method retrieves inflated LR") +
    facet_wrap(~ FC_nReceiverCells)
  
  ##########################################
  ##### Generate plots for L_parameter #####
  
  
  
  ##########################################################
  ##### generate precision/recall plots across indexLR #####
  
  # the goal here is to understand how values of recall and precision change depending of ligand/receptor pair (indexLR) that one chooses
  # We are calculating how many cells express ligand/receptor in CT1/CT2 compared to all cells expressing ligand/receptor regardless of celltype
  
  # For the final plot, too many combinations of parameters are used, reduce them here
  n = length(FC_nSenderCells)
  FC_nSenderCells_touse = FC_nSenderCells[unique(c(seq(1, n, by = 2), n))]
  n = length(FC_nReceiverCells)
  FC_nReceiverCells_touse = FC_nReceiverCells[unique(c(seq(1, n, by = 2), n))]
  
  # here we are summarising the data where we merge information of how many cells per celltype and combination of parameters are expressing ligands/receptors 
  # and the f1score/precision/recall scores for further plotting
  lst_data = list()
  for(indexLR in config$indexLR_toSample)
  {
    # generate a dataframe of cell number statistics for CT1 (expressing ligands)
    ligands = lapply(gene_metadata_lst[[dataset]] , function(x) {x[["CT1"]] %>% do.call(rbind.data.frame,.)
    }) %>%  
      do.call(rbind.data.frame,.) %>%
      tibble::rownames_to_column("filename") %>%
      mutate(filename = gsub("(_indexLR_).*","",filename)) %>%
      filter(celltype == "CT1") %>%
      dplyr::group_by(indexLR,filename) %>%
      dplyr::summarise(
        across(where(is.numeric), ~mean(.x, na.rm = TRUE)),           
        across(where(is.character), ~paste(.x, collapse = ", "))      
      ) %>% mutate(filename = str_c(filename, "_indexLR_",indexLR))
    
    # generate a dataframe of cell number statistics for CT2 (expressing receptors)
    receptors = lapply(gene_metadata_lst[[dataset]] , function(x) {x[["CT2"]] %>% do.call(rbind.data.frame,.)
    }) %>% 
      do.call(rbind.data.frame,.) %>%
      tibble::rownames_to_column("filename") %>%
      mutate(filename = gsub("(_indexLR_).*","",filename)) %>%
      filter(celltype == "CT2") %>%
      dplyr::group_by(indexLR,filename) %>%
      dplyr::summarise(
        across(where(is.numeric), ~mean(.x, na.rm = TRUE)),           # Numeric: Mean
        across(where(is.character), ~paste(.x, collapse = ", "))      # String: Concatenate
      ) %>% mutate(filename = str_c(filename, "_indexLR_",indexLR))
    
    
    # iterate over every method and add f1score/precision/recall values for each combination of FC parameters
    for(m in names(statistics_results_lst_recallprecision_non_indexLRaveraged[[dataset]]))
    {
      tmp_df = statistics_results_lst_recallprecision_non_indexLRaveraged[[dataset]][[m]] %>%
        select(precision,recall,FC_nSenderCells,FC_nReceiverCells,filename) %>%
        filter(FC_nSenderCells %in% FC_nSenderCells_touse & FC_nReceiverCells %in% FC_nReceiverCells_touse) %>%
        mutate(filename = str_remove(filename, "_l_\\d+")) %>%
        dplyr::group_by(filename) %>%
        dplyr::summarise(
          across(where(is.numeric) | is.logical, ~mean(.x, na.rm = TRUE))
        )
      
      tmp_df2 = left_join(ligands, receptors, by = "filename")
      tmp_df3 = left_join(tmp_df,tmp_df2 ,by = "filename")
      colnames(tmp_df3) %<>% gsub("[.]x", "_ligand", .) %>% gsub("[.]y", "_receptor", .)
      
      lst_data[[m]] = tmp_df3
    }
  }
  
  # summarise data based on indexLR to calculate mean and standard deviation
  summary_table = lapply(lst_data, function(x) {
    x %>% mutate(filename_without_indexLR = str_remove(filename, "_indexLR_\\d+"), indexLR = as.factor(x$indexLR_ligand)) # ligand and receptor indexLR is the same 
  })
  
  # combine all methods into a single dataframe
  final_df = bind_rows(summary_table, .id = "method")
  final_df$filename_without_indexLR %<>% gsub("FC_nSenderCells_", "S" ,.) %>% gsub("FC_nReceiverCells_", "R",.)
  
  # Use regex to simplify the filename column for easier plotting
  final_df$filename_without_indexLR <- str_replace(final_df$filename_without_indexLR, "(?<=_R)(\\d+\\.?\\d*)", function(m) {
    sprintf("%.1f", as.numeric(m))
  }) %>% str_replace(., "(?<=S)(\\d+\\.?\\d*)", function(m) {
    sprintf("%.1f", as.numeric(m))
  }) %>% str_replace_all(., "\\.0", "")
  
  # create factors
  final_df$filename_without_indexLR = factor(final_df$filename_without_indexLR, levels = c("S2_R2", "S2_R6", "S2_R10","S6_R2", "S6_R6","S6_R10", "S10_R2", "S10_R6","S10_R10"))
  
  indexLR_plot1 = ggplot(final_df, aes(x = precision, y = filename_without_indexLR, fill = filename_without_indexLR)) +
    geom_density_ridges(stat = "binline", bins = 30, scale = 0.9, alpha = 0.85, draw_baseline = FALSE) +
    facet_wrap(~method, ncol = 2) +
    scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
    theme_bw() + 
    labs(
      title = "Method performance across different indexLR",
      x = "Precision",
      y = "PC Snd/Rcv parameters"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  indexLR_plot2 = ggplot(final_df, aes(x = recall, y = filename_without_indexLR, fill = filename_without_indexLR)) +
    geom_density_ridges(stat = "binline", bins = 30, scale = 0.9, alpha = 0.85, draw_baseline = FALSE) +
    facet_wrap(~method, ncol = 2) +
    scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
    theme_bw() + 
    labs(
      title = "Method performance across different indexLR",
      x = "Recall",
      y = "PC Snd/Rcv parameters"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  ######################################################
  ##### save precision/recall plots across indexLR #####
  ######################################################
  
  pdf(file.path(path_results_dir ,paste0(dataset, "_recall_precision_plots.pdf")), width = 15, height = 7)
  p1 %>% print
  ggarrange(plotlist = list(p3,p4), common.legend = T) %>% print
  heatmap1 %>% print
  heatmap2_faceted_f1score %>% print
  heatmap2_faceted_recall %>% print
  ridge_plot %>% print
  p2 %>% print
  N_interactions_f1score_plot %>% print
  upset1 %>% print
  upset2 %>% print
  indexLR_plot1 %>% print
  indexLR_plot2 %>% print
  dev.off()
  
  
  pdf(file.path(path_results_dir ,paste0(dataset, "_recall_precision_plots2.pdf")), width = 20, height = 7)
  heatmap2 %>% print
  heatmap3 %>% print
  heatmap4 %>% print
  dev.off()
}


# ==========================================
# Generate a summary plot across datasets
# ==========================================
tmp_df_long = lapply(names(statistics_results_lst_recallprecision), function(dataset) {
  do.call(rbind.data.frame, statistics_results_lst_recallprecision[[dataset]]) %>% 
    mutate(dataset = dataset)
}) %>% 
  do.call(rbind.data.frame, .) %>%
  mutate(
    FC_nReceiverCells = as.numeric(as.character(FC_nReceiverCells)),
    FC_nSenderCells = as.numeric(as.character(FC_nSenderCells))
  ) %>%
  arrange(FC_nReceiverCells, FC_nSenderCells) %>%
  mutate(params = paste0("FC_nReceiverCells_", FC_nReceiverCells, "_FC_nSenderCells_", FC_nSenderCells, "_l_",l_param_index)) %>%
  mutate(params = factor(params, levels = unique(params)))

'
df_relative <- tmp_df_long %>%
  group_by(dataset, params) %>%
  dplyr::mutate(relative_score = recall / sum(recall, na.rm = TRUE)) %>%
  dplyr::summarise(
    rel_vec = list(setNames(relative_score, method)),
    .groups = "drop"
  )

mat_rel <- df_relative %>%
  pivot_wider(names_from = params, values_from = rel_vec) %>%
  tibble::column_to_rownames("dataset") %>%
  as.matrix()
'

plot_data <- tmp_df_long %>%
  mutate(
    # Step 1: Extract the numerical digits directly out of the parameter string
    rcv_val = str_split_i(params, "_", 3) %>% as.numeric(),
    snd_val = str_split_i(params, "_", 6) %>% as.numeric(),
    
    # Step 2: Calculate the product feature for the continuous X-axis
    Snd_Rcv_Product = rcv_val * snd_val,
    
    # Step 3: Format faceting variables as clean factors
    dataset = factor(dataset),
    l_param_index = factor(paste0("l_param:", l_param_index)),
    
    # Optional: Build your Boolean metadata flag for spatial awareness
    AccountsSpatialDistance = !method %in% c("cellphonedb", "seurat_wilcoxon")
  )

# Define manual dashed patterns for your baseline metrics
all_methods <- unique(plot_data$method)
line_types  <- setNames(rep("solid", length(all_methods)), all_methods)
line_types["cellphonedb"]    <- "dashed"
line_types["seurat_wilcoxon"] <- "dashed"

product_f1score_plot = ggplot(plot_data, aes(x = Snd_Rcv_Product, y = f1score, 
                                             color = method, 
                                             linetype = method,
                                             linewidth = AccountsSpatialDistance,
                                             group = method)) +
  
  geom_smooth(method = "loess", 
              span = 0.75, 
              se = FALSE, 
              alpha = 0.85) +
  
  facet_grid(dataset ~ l_param_index, scales = "free_x") + 
  
  scale_color_viridis_d(option = "D", name = "Methods") +
  scale_linetype_manual(values = line_types, name = "Methods") +
  
  # Keep the uniform linewidths as you configured
  scale_linewidth_manual(values = c("TRUE" = 0.5, "FALSE" = 0.5), 
                         name = "AccountsSpatialDistance") +
  
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  
  labs(
    x = "Combined parameter density (Sender × Receiver Product Scale)",
    y = "F1score Value",
    title = "Smoothed f1score trends across parameter"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray95"),
    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.4, "lines"),
    legend.key.width = unit(1.5, "cm"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 9)
  ) + 
  guides(
    linetype = guide_legend(override.aes = list(linewidth = 0.5)),
    linewidth = guide_legend(override.aes = list(
      linetype = c("TRUE" = "dashed", "FALSE" = "solid"),
      color = "black" # Force key icons to be readable black lines instead of missing colors
    ))
  )

ggsave(filename = "output/results/summary_across_datasets.png", plot = product_f1score_plot, width = 200, height = 150, units = "mm")


'
# make plot for precision/recall for colocalization and original
out = lapply(datasets, function(dataset) {
  message(dataset)
  # NON-COLOCALIZATION
  metric_results = readRDS("output_recent/final_scores.RDS")
  statistics_results_lst_recallprecision_non_colocalization = metric_results %>% lapply(., function(dataset) {
    # iterate over each method and average metrics of indexLR parameter
    lapply(dataset, function(method) {
      tmp_df = method %>% do.call(rbind, .) %>%
        # Create the new column by removing the suffix
        # The regex "_indexLR_\\d+$" targets "_indexLR_" and all digits at the end
        mutate(file_without_indexLR = str_remove(rownames(.), "_indexLR_\\d+$")) %>%
        group_by(file_without_indexLR) %>%
        # average every numeric column
        dplyr::summarise(
          across(
            everything(), 
            ~ if(is.numeric(.x)) mean(.x, na.rm = TRUE) else .x[1]
          ), 
          .groups = "drop"
        )
    })
  })
  
  # COLOCALIZATION
  metric_results = readRDS("output_colocalization/final_scores.RDS")
  statistics_results_lst_recallprecision_colocalization = metric_results %>% lapply(., function(dataset) {
    # iterate over each method and average metrics of indexLR parameter
    lapply(dataset, function(method) {
      tmp_df = method %>% do.call(rbind, .) %>%
        # Create the new column by removing the suffix
        # The regex "_indexLR_\\d+$" targets "_indexLR_" and all digits at the end
        mutate(file_without_indexLR = str_remove(rownames(.), "_indexLR_\\d+$")) %>%
        group_by(file_without_indexLR) %>%
        # average every numeric column
        dplyr::summarise(
          across(
            everything(), 
            ~ if(is.numeric(.x)) mean(.x, na.rm = TRUE) else .x[1]
          ), 
          .groups = "drop"
        )
    })
  })
  
  data_across_methods_dispersedSignal = statistics_results_lst_recallprecision_non_colocalization[[dataset]] %>% do.call(rbind, .)
  data_across_methods_colocalization = statistics_results_lst_recallprecision_colocalization[[dataset]] %>% do.call(rbind, .)
  # Assuming both dataframes have columns: "method", "recall", and "precision"
  data_across_methods_dispersedSignal$condition <- "Dispersed signal"
  data_across_methods_colocalization$condition <- "Colocalization"
  data_across_methods_dispersedSignal$dataset <- dataset
  data_across_methods_colocalization$dataset <- dataset
  
  combined_data <- bind_rows(data_across_methods_dispersedSignal, data_across_methods_colocalization)
  
  long_data <- combined_data %>%
    pivot_longer(cols = c(recall, precision), # Change to precision or false_positives based on your columns
                 names_to = "Metric", 
                 values_to = "Value") %>%
    mutate(Metric = recode(Metric, 
                           "recall" = "Recall", 
                           "precision" = "Precision"))
  
  return(long_data)  
}) %>% setNames(datasets)


long_data_plot = do.call(rbind, out)

# 1. Pool data across your 3 datasets carefully preserving "Metric"
pooled_data <- long_data_plot %>%
  # Force Metric to be exactly what ggplot expects
  dplyr::group_by(method, condition, Metric, FC_nSenderCells, FC_nReceiverCells, l_param_index) %>%
  dplyr::summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# 2. Run the simplified plot
p = ggplot(pooled_data, aes(x = method, y = Value, color = condition)) +
  stat_summary(
    fun.data = median_hilow, 
    fun.args = list(conf.int = 0.5),
    geom = "pointrange",
    position = position_dodge(width = 0.5),
    size = 0.7,
    linewidth = 0.9
  ) +
  
  # This matches the capital "Metric" column from above
  facet_wrap(~Metric, scales = "free_y", nrow = 2) + 
  
  scale_color_manual(values = c("#4292C6", "#EF3B2C")) + 
  labs(
    title = "Benchmark performance profile: Dispersed signal vs Colocalization",
    subtitle = "Aggregated across spatial datasets (CosMx, MERFISH, Visium_HD)",
    x = "Inference Method",
    y = "Score",
    color = "Condition"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(filename = "output/results/summary_across_datasets_dispersedSignal_vs_colocalization.png", plot = p, width = 200, height = 150, units = "mm")
'
