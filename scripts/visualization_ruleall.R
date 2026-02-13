suppressMessages({
  # Load package
  library(dplyr)
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
  library(ggspavis)
  library(sf)
  library(fmsb)
  source("scripts/helper_functions.R")
})

# Get list with command line arguments by name
option_list = list(
  make_option(c("--config.yaml"), type="character", default=NULL, help="config file", metavar="character"),
  make_option(c("--path_output_dir"), type="character", default=NULL, help="path to output directory", metavar="character"),
  make_option(c("--path_results_dir"), type="character", default=NULL, help="path to where to save the results", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# An useful variability if the argument is missing
if (is.null(opt$config.yaml) | is.null(opt$path_output_dir) | is.null(opt$path_results_dir)){
  print_help(opt_parser)
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Call the argument
path_config.yaml = opt$config.yaml
path_output_dir = opt$path_output_dir
path_results_dir = opt$path_results_dir
dir.create(path_results_dir)

########################################
##### Loading and processing files #####
########################################

metric_results = readRDS("output/final_scores.RDS")

config = yaml::read_yaml(path_config.yaml)

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
diagnostic_df_realFC = list()
diagnostic_plots_perCT = list()
cumulative_expression_lst = list()
for(dataset in datasets)
{
  
  # Load visium data to generate a plot
  if(FALSE & dataset %in% c("Visium_HPC_SJ"))
  {
    spe = read10xVisium(file.path("data", dataset)) 
    rownames(spe) = rowData(spe)$symbol %>% toupper()
    colData(spe) %<>% cbind(., master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$metadata)
  }
  
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
    inflated_counts_file = inflated_counts_files %>% extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_indexLR_",x["indexLR"] , ".tsv$"), inflated_counts_files))
    inflated_counts = read.table(file.path(file_path,inflated_counts_file)) 
    
    # load correct simulated_cellmetadata file depending on the params_grid
    simulated_cellmetadata_file = simulated_cellmetadata_files %>% extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_indexLR_",x["indexLR"] , ".json$"), simulated_cellmetadata_files))
    master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]] = read_json(file.path(file_path,simulated_cellmetadata_file)) %>% convert_json_to_df
    
    # load correct simulated interactions file depending on the params_grid
    simulated_interactions_file = simulated_interactions_files %>% extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_indexLR_",x["indexLR"] , ".RDS$"), simulated_interactions_files))
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
        
        cumulative_expression_lst[[naming]][[CT]][[gene]] = data.frame(
          expression = inflated_counts[gene, ] %>% as.vector %>% unlist,
          celltype = m$Celltype 
        ) %>%
          group_by(celltype) %>%
          dplyr::summarise(
            # A cell expresses the gene if count > 0
            num_expressing = sum(expression > 0),
            total_cells = n(),
            pct_expressing = round(100 * (num_expressing / total_cells), 2),
            indexLR = x["indexLR"]
          )
        
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
  ############################################################################
  ##### Plot cumulative gene expression by radius for all receptor genes #####
  LR_database = read.table("data/LR_database.tsv")
  genes = LR_database$receptor %>% str_split("_") %>% unlist %>% unique
  genes = genes[genes %in% rownames(inflated_counts)]
  
  cumulative_expression_acrossRadius = list()
  lst_pvals = list()
  lst = list()
  for(gene in genes)
  {
    # split the data into groups of 50 cells for each celltype to compute statistics in the end
    seeds = seq(1:20)
    for(seed in seeds)
    {
      x = cumulative_expression(gene = gene,seed = seed, counts = original_counts, cellmetadata = master_lst_diagnosticPlots[[1]]$simulated_cellmetadata, size = 50)
      cumulative_expression_acrossRadius[[gene]][[paste0("seed_",seed)]] = x$plot
      lst[[gene]][[paste0("seed_",seed)]] = x$distance_curve_ct %>% mutate(seed = paste0("seed_",seed))
    }
    # compute statistics for each radius
    df = do.call(rbind.data.frame, lst[[gene]])
    
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
}

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
  
  barplot(N_genes_CT1CT2_significant$n_significant,
          names.arg = N_genes_CT1CT2_significant$radius,
          xlab = "radius",
          ylab = "count of significant results",
          main = "Amount of Receptors for CT1-CT2 whose cumulative expression is statistically significant for each radius",
          col = "lightblue") %>% print
  
  ggarrange(plotlist = list(cumulative_expression_acrossRadius[[1]]$seed_1, 
                            cumulative_expression_acrossRadius[[2]]$seed_5 , 
                            cumulative_expression_acrossRadius[[3]]$seed_8,
                            cumulative_expression_acrossRadius[[4]]$seed_15)) %>% print
  dev.off()
}

##################################
##### precision/recall plots #####
##################################

# plot TPR/sensitivity/recall
statistics_results_lst_recallprecision_plot = list()
master_lst_precision_recall = list()

for(dataset in datasets)
{
  for(method in methods)
  {
    # Load files for upset plot
    significant_interactions_files = file.path(path_output_dir,dataset,method) %>% 
      list.files(., pattern = "significant_interactions")
    
    for(file in significant_interactions_files)
    {
      master_lst_precision_recall[["significant_interactions"]][[method]][[file]] = read.table(file.path(path_output_dir,dataset,method,file),header = T) %>% 
        mutate(significant = as.logical(significant)) %>%
        filter(significant) %>% 
        select(ligand_receptor) %>% 
        unlist %>%
        unname
    }
    
    
    ################################################################################################
    ##### Average the metrics for different indexLR but across same combination of parameters  #####
    
    statistics_results_lst_recallprecision_plot_non_indexLRaveraged = metric_results %>% lapply(., function(dataset) {
      # iterate over each method and average metrics of indexLR parameter
      lapply(dataset, function(method) {
        tmp_df = method %>% do.call(rbind, .) %>%
          # Create the new column by removing the suffix
          # The regex "_indexLR_\\d+$" targets "_indexLR_" and all digits at the end
          mutate(filename = rownames(.))
      })
    })
    
    # extract rows that have the same parameters besides indexLR and then average the columns. We are averaging results of indexLR
    # generate the grid of parameters used for naming the list to generate output
    # this params_grid is different from previous because before we didnt look at methods and now we have different l_parameters for each method
    statistics_results_lst_recallprecision_plot = metric_results %>% lapply(., function(dataset) {
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
  
  ###########################################################################
  ##### Calculate % of cells from cellmetadata that we added signal to  #####
  # we dont care about l_param_index parameter (here we use it for matching with master list) nor indexLR as that doesnt affect amount of cells meaning all files in simulated_cellmetadata_file should contain same %
  # load correct simulated_cellmetadata file depending on the params_grid
  
  params_grid = expand.grid(FC_nSenderCells = FC_nSenderCells, FC_nReceiverCells = FC_nReceiverCells, l_param_index = l_param_index)
  tmp_lst = list()
  for(row in 1:nrow(params_grid))
  {
    x = params_grid[row,] %>% as.numeric ;  names(x) = c("FC_nSenderCells","FC_nReceiverCells", "l_param_index")
    
    naming = paste0("FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_l_", x["l_param_index"])
    
    simulated_cellmetadata_file = simulated_cellmetadata_files %>% extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"]), simulated_cellmetadata_files))
    
    y = read_json(file.path(file_path,simulated_cellmetadata_file[1])) %>% convert_json_to_df # take only 1 file -> to speed up computations
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
  reordered_list = tmp_lst[statistics_results_lst_recallprecision_plot[[dataset]][[method]]$file_without_indexLR]
  
  # make sure the entire reordered vector has same names as the dataframe column
  reordered_list_check = reordered_list %>% names %>% is.na
  stopifnot(!any(reordered_list_check))
  
  #### Combine all methods together
  # for precision recall plots
  # add pct_signal_added variable
  statistics_results_lst_recallprecision_plot[[dataset]] = do.call(rbind, statistics_results_lst_recallprecision_plot[[dataset]]) %>%
    mutate(pct_CT1cells_expressing_L = 
             rep(reordered_list %>% do.call(rbind.data.frame,.) %>% filter(Celltype == "CT1") %>% pull(pct_signal_added), length(methods)),
           N_CT1cells_expressing_L = 
             rep(reordered_list %>% do.call(rbind.data.frame,.) %>% filter(Celltype == "CT1") %>% pull(amount_signal_added), length(methods))
    )
  statistics_results_lst_recallprecision_plot[[dataset]] = statistics_results_lst_recallprecision_plot[[dataset]] %>%
    mutate(pct_CT2cells_expressing_R = 
             rep(reordered_list %>% do.call(rbind.data.frame,.) %>% filter(Celltype == "CT2") %>% pull(pct_signal_added), length(methods)),
           N_CT2cells_expressing_R = 
             rep(reordered_list %>% do.call(rbind.data.frame,.) %>% filter(Celltype == "CT2") %>% pull(amount_signal_added), length(methods))
    )
  
  ###########################
  ##### Plot showing FP #####
  
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  FP = ggplot(tmp_df, aes(y = FP, x = method, color = method)) +
    geom_point() +
    scale_y_log10() +
    ylab("FP") +
    ggtitle("Number of FP per method across all combination of parameters")
  
  
  #######################
  ##### UpSet plots #####
  
  interactions_retrieved = lapply(master_lst_precision_recall[["significant_interactions"]] , function(method) {
    # find all unique elements
    all_strings = unique(unlist(master_lst_precision_recall$significant_interactions))
    # For each element in the list, check if the global strings exist there
    binary_df = as.data.frame(lapply(method, function(x) {
      as.numeric(all_strings %in% x)
    }))
    
    data.frame(method = ifelse(rowSums(binary_df) > 0, 1,0), row.names = all_strings)
  }) %>% do.call(cbind,.)
  
  colnames(interactions_retrieved) = names(master_lst_precision_recall$significant_interactions)
  
  upset = upset(
    interactions_retrieved, 
    colnames(interactions_retrieved),
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        mapping = aes(fill = "bars") # You can style the bars here
      )
    ),
    set_sizes = (
      upset_set_size() + 
        # Use expand_limits to ensure the axis goes high enough for the labels
        expand_limits(y = 1500) +
        geom_text(aes(label = ..count..), hjust = 1.1, stat = 'count') +
        expand_limits(y = 120)
    ),
    themes = upset_default_themes(text = element_text(size = 12),
                                  legend.position = "none")
  )
  
  ##########################################
  ##### Generate Precision recall plot #####
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]] %>%
    mutate(ratio_ReceiverSender = (FC_nReceiverCells/FC_nSenderCells) %>% log2) %>% 
    group_by(ratio_ReceiverSender,method ) %>%
    summarise_at(vars(precision,recall), mean)
  
  # dont average across ratio_ReceiverSender
  tmp_df2 = statistics_results_lst_recallprecision_plot[[dataset]] %>%
    mutate(ratio_ReceiverSender = (FC_nReceiverCells/FC_nSenderCells) %>% log2)
  
  p1 = ggplot(tmp_df2, aes(x = precision, y = recall)) +
    geom_density_2d(
      aes(fill = after_stat(level)), 
      contour_var = "ndensity",
      h = c(0.1,0.1) # fine tude KDE bandwidth as some methods have very low variance and density estimation doesnt work
    ) +
    facet_wrap(~ method, ncol = 4) +
    labs(title = "Independent 2D Density per method")
  
  ###################################
  ##### Generate f1 score plots #####
  
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  
  
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
  
  p3 = ggplot(tmp_df, aes(x = pct_CT1cells_expressing_L, y = f1score, color = method)) + 
    geom_point(size = 0.5) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l_param_index)+
    ggtitle("F1score faceted by L") 
  
  p4 = ggplot(tmp_df, aes(x = N_CT1cells_expressing_L, y = f1score, color = method)) + 
    geom_point(size = 0.5) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l_param_index)+
    ggtitle("F1score faceted by L") 
  
  p5 = ggplot(tmp_df, aes(x = pct_CT2cells_expressing_R, y = f1score, color = method)) + 
    geom_point(size = 0.5) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l_param_index)+
    ggtitle("F1score faceted by L") 
  
  p6 = ggplot(tmp_df, aes(x = N_CT2cells_expressing_R, y = f1score, color = method)) + 
    geom_point(size = 0.5) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l_param_index)+
    ggtitle("F1score L") 
  
  ##########################################
  ##### Generate ranking LR genes plot #####
  
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]] %>% 
    mutate(rank = round(log10(rank),2))
  
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
  
  tmp_df2 = tmp_df %>%
    na.omit %>%
    group_by(method,l_param_index) %>%
    summarise_at(vars(rank), list(rank = mean)) %>%
    mutate(rank = round(rank,2))
  
  # Heatmap of rank according to l_param_index radius
  heatmap2 = ggplot(tmp_df2, aes(method, l_param_index, fill= rank)) +
    geom_tile(color = "black") +
    geom_text(aes(label = rank), color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("rank of how each method retrieves inflated LR")+
    scale_y_continuous(breaks=l_param_index) +
    ylab("Index of l_param_index parameter") + 
    labs(fill = "log10(rank)")
}

##########################################################
##### generate precision/recall plots across indexLR #####
##########################################################
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
  ligands = lapply(cumulative_expression_lst , function(x) {x[["CT1"]] %>% do.call(rbind.data.frame,.) %>%
      mutate(perc_expressing_ofTotalexpressed = num_expressing/sum(num_expressing)) 
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
  receptors = lapply(cumulative_expression_lst , function(x) {x[["CT2"]] %>% do.call(rbind.data.frame,.) %>%
      mutate(perc_expressing_ofTotalexpressed = num_expressing/sum(num_expressing)) 
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
  for(m in names(statistics_results_lst_recallprecision_plot_non_indexLRaveraged$Visium_HD_HPC))
  {
    tmp_df = statistics_results_lst_recallprecision_plot_non_indexLRaveraged$Visium_HD_HPC[[m]] %>%
      filter(FC_nSenderCells %in% FC_nSenderCells_touse & FC_nReceiverCells %in% FC_nReceiverCells_touse) %>%
      mutate(filename = str_remove(filename, "_l_\\d+")) %>%
      dplyr::group_by(filename) %>%
      dplyr::summarise(
        across(where(is.numeric), ~mean(.x, na.rm = TRUE)),           # Numeric: Mean
        across(where(is.character), ~paste(.x, collapse = ", "))      # String: Concatenate
      )
    
    tmp_df2 = left_join(ligands, receptors, by = "filename")
    tmp_df3 = left_join(tmp_df,tmp_df2 ,by = "filename")
    colnames(tmp_df3) %<>% gsub("[.]x", "_ligand", .) %>% gsub("[.]y", "_receptor", .)

    lst_data[[m]] = tmp_df3
  }
}

# summarise data based on indexLR to calculate mean and standard deviation
summary_table = lapply(lst_data, function(x) {
  x %>% 
    mutate(filename_without_indexLR = str_remove(filename, "_indexLR_\\d+"), indexLR <- as.factor(x$indexLR)) %>%
    group_by(filename_without_indexLR) %>%
    dplyr::summarise(
      mean_precision = mean(precision),
      sd_precision = sd(precision),
      mean_recall = mean(recall),
      sd_recall = sd(recall),
      perc_expressing_ofTotalexpressed_ligand_overall = mean(perc_expressing_ofTotalexpressed_ligand) %>% round(3),
      perc_expressing_ofTotalexpressed_receptor_overall = mean(perc_expressing_ofTotalexpressed_receptor %>% round(3))
    ) %>% ungroup
})

# combine all methods into a single dataframe
final_df = bind_rows(summary_table, .id = "method")
final_df$filename_without_indexLR %<>% gsub("FC_nSenderCells_", "S" ,.) %>% gsub("FC_nReceiverCells_", "R",.) 
final_df %<>% mutate(filename_with_perc_cells_expressing = str_c("Total_ligand_", round(perc_expressing_ofTotalexpressed_ligand_overall,3),"_",
                                                                 "Total_receptor_", round(perc_expressing_ofTotalexpressed_receptor_overall,3)))

# Use regex to simplify the filename column for easier plotting
final_df$filename_without_indexLR <- str_replace(final_df$filename_without_indexLR, "(?<=_R)(\\d+\\.?\\d*)", function(m) {
  sprintf("%.1f", as.numeric(m))
}) %>% str_replace(., "(?<=S)(\\d+\\.?\\d*)", function(m) {
  sprintf("%.1f", as.numeric(m))
})

indexLR_plot1 = ggplot(final_df, aes(x = filename_without_indexLR, y = mean_precision, color = method)) +
  geom_errorbar(aes(ymin = mean_precision - sd_precision, 
                    ymax = mean_precision + sd_precision), 
                width = 0.2) +
  geom_point(size = 3) +
  coord_flip() + # Flip the coordinates to make filenames readable
  facet_wrap(~method) +
  theme_bw() +
  labs(
    title = "Precision by FC param combination",
    subtitle = "Error bars represent ±1 Standard Deviation",
    x = "Filename",
    y = "Averaged across indexLR precision value"
  )

indexLR_plot2 = ggplot(final_df, aes(x = filename_without_indexLR, y = mean_recall, color = method)) +
  geom_errorbar(aes(ymin = mean_recall - sd_recall, 
                    ymax = mean_recall + sd_recall), 
                width = 0.2) +
  geom_point(size = 3) +
  coord_flip() +
  facet_wrap(~method) +
  theme_bw() +
  labs(
    title = "Recall by FC param combination",
    subtitle = "Error bars represent ±1 Standard Deviation",
    x = "Filename",
    y = "Averaged across indexLR recall value"
  )

######################################################
##### save precision/recall plots across indexLR #####
######################################################
pdf(file.path(path_results_dir ,paste0(dataset, "_recall_precision_plots.pdf")), width = 12, height = 7)
p1 %>% print
ggarrange(plotlist = list(p3,p4,p5,p6), common.legend = T) %>% print
heatmap1 %>% print
heatmap2 %>% print
p2 %>% print
FP %>% print
upset %>% print
indexLR_plot1 %>% print
indexLR_plot2 %>% print
# Generate a clean table describing FC params to real values of how many cells out of all cells in the dataset express ligand/receptor
knitr::kable(head(final_df) %>% select(filename_without_indexLR, perc_expressing_ofTotalexpressed_ligand_overall,perc_expressing_ofTotalexpressed_receptor_overall)
             , caption = "FC params to % of all cells in the dataset that express L/R ; averaged across all L/R genes")
dev.off()