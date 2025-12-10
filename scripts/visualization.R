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
for(dataset in datasets)
{
  
  # LOad visium data to generate a plot
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
    
    '
    ###########################################################
    ##### Generate Plots of real FC after semi-simulation #####
    
    # Select file realFC after simulation
    #realFC_aftersimulation_file = file_path %>% list.files(., pattern = "realFC_aftersimulation") %>% extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_indexLR_",x["indexLR"], ".RDS$"), .))
    realFC_aftersimulation_file = file_path %>% list.files(., pattern = "realFC_aftersimulation") %>% extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_indexLR","_FC_nReceiverCells_", x["FC_nReceiverCells"],"_",x["indexLR"], ".RDS$"), .))
    
    tmp_lst = readRDS(file.path(file_path,realFC_aftersimulation_file)) %>% 
      unlist %>% subset(.,!is.infinite(.)) %>% 
      plot_FCafter_semisimulation(. , indexLR = x["indexLR"] ,theoreticalFC = FC, FC_nSenderCells = x["FC_nSenderCells"], FC_nReceiverCells = x["FC_nReceiverCells"])
    
    diagnostic_plots_realFC[[dataset]][[paste0("indexLR_",x["indexLR"])]][[paste0("n_neighbors_",x["n_neighbors"])]] = tmp_lst %>% pluck("plot")
    
    diagnostic_df_realFC[[dataset]][[paste0("index_",i)]] = tmp_lst[c(2,3,4)] %>% as.data.frame
    
    '
    
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
      x = cumulative_expression(gene = gene,seed = seed, counts = original_counts, cellmetadata = master_lst_diagnosticPlots[[1]]$simulated_cellmetadata)
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
                            cumulative_expression_acrossRadius[[4]]$seed_15, 
                            cumulative_expression_acrossRadius[[5]]$seed_12,
                            cumulative_expression_acrossRadius[[6]]$seed_19))
  dev.off()
}

##################################
##### precision/recall plots #####
##################################

# plot TPR/sensitivity/recall
statistics_results_lst_recallprecision_plot = list()
master_lst_precision_recall = list()
ranking_LRgenes_lst_plot = list()
for(dataset in datasets)
{
  for(method in methods)
  {
    # load metric files
    metrics_files = file.path(path_output_dir,dataset,"metrics/",method) %>%
      list.files(., pattern = "f1_score")
    ranking_LRgenes_files = file.path(path_output_dir,dataset,"metrics/",method) %>% 
      list.files(., pattern = "ranking_LRgenes")
    
    for(file in metrics_files)
    {
      master_lst_precision_recall[["metrics"]][[file]] = read.table((file.path(path_output_dir,dataset,"metrics/",method,file)), header = TRUE)
    }
    
    # NAs here mean that the LR gene that we inflated is not present in the output of the method
    for(file in ranking_LRgenes_files)
    {
      master_lst_precision_recall[["ranking"]][[file]] = read.table((file.path(path_output_dir,dataset,"metrics/",method,file))) %>% unlist %>% as.numeric
    }
    
    
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] = do.call(rbind,master_lst_precision_recall[["metrics"]]) %>% as.data.frame
    ranking_LRgenes_lst_plot[[dataset]][[method]] = do.call(rbind,master_lst_precision_recall[["ranking"]]) %>% as.data.frame
    
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
    
    # extract rows that have the same parameters besides indexLR and then average the columns. We are averaging results of indexLR
    # generate the grid of parameters used for naming the list to generate output
    # this params_grid is different from previous because before we didnt look at methods and now we have different l_parameters for each method
    params_grid = expand.grid(FC_nSenderCells = FC_nSenderCells, FC_nReceiverCells = FC_nReceiverCells, l = l_param_index)
    
    tmp_lst = list()
    for(row in 1:nrow(params_grid))
    {
      x = params_grid[row,] %>% as.numeric ;  names(x) = c("FC_nSenderCells","FC_nReceiverCells", "l")
      
      naming = paste0("FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_l_", x["l"])
      
      # extract filenames -> compute mean across columns and save to temporary list
      filenames = statistics_results_lst_recallprecision_plot[[dataset]][[method]] %>% 
        rownames() %>% 
        extract(grepl(naming, .))
      
      tmp_lst[[row]] = statistics_results_lst_recallprecision_plot[[dataset]][[method]][filenames,] %>% colMeans()
      
      ###########################################################################
      ##### Calculate % of cells from cellmetadata that we added signal to  #####
      
      # we dont care about l parameter nor indexLR as that doesnt affect amount of cells meaning all files in simulated_cellmetadata_file should contain same %
      # load correct simulated_cellmetadata file depending on the params_grid
      simulated_cellmetadata_file = simulated_cellmetadata_files %>% extract(grepl(paste0("_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"]), simulated_cellmetadata_files))
      
      y = read_json(file.path(file_path,simulated_cellmetadata_file[1])) %>% convert_json_to_df # take only 1 file -> to speed up computations
      y = y$metadata
      # calculate percentage of "CT1_signal_added" and "CT2_signal_added" cells
      pct_signal_added = y %>% mutate(
        parent = case_when(
          Celltype_updated %in% c("CT1", "CT1_signalAdded") ~ "CT1",
          Celltype_updated %in% c("CT2", "CT2_signalAdded") ~ "CT2"
        )
      ) %>%
        filter(!is.na(parent)) %>% 
        group_by(parent) %>%
        summarise(
          pct_signal_added = round(100 * mean(grepl("signalAdded", Celltype_updated)))
        ) 
      
      # compute mean across indexLR and add % cellsl
      tmp_lst[[row]] = c(tmp_lst[[row]] , pct_signal_added) %>% unlist
      names(tmp_lst)[row] = naming
    }
    
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] = do.call(rbind,tmp_lst) %>% as.data.frame
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] %<>% mutate(FC_nSenderCells = params_grid$FC_nSenderCells, 
                                                                                 FC_nReceiverCells = params_grid$FC_nReceiverCells,
                                                                                 l = params_grid$l,
                                                                                 method = method,
                                                                                 dataset = dataset)
    ###############################
    ##### Ranking of LR genes #####
    
    # extract rows that have the same parameters besides indexLR and then average the columns. We are averaging results of indexLR
    # generate the grid of parameters used for naming the list to generate output
    tmp_lst = list()
    
    for(row in 1:nrow(params_grid))
    {
      x = params_grid[row,] %>% as.numeric ;  names(x) = c("FC_nSenderCells","FC_nReceiverCells", "l")
      
      naming = paste0("FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_l_", x["l"])
      
      # extract filenames -> compute mean across columns and save to temporary list
      filenames = ranking_LRgenes_lst_plot[[dataset]][[method]] %>% 
        rownames() %>% 
        extract(grepl(naming, .))
      
      tmp_lst[[row]] = ranking_LRgenes_lst_plot[[dataset]][[method]][filenames,] %>% mean(., na.rm = TRUE)
      
      names(tmp_lst)[row] = naming
    }
    
    ranking_LRgenes_lst_plot[[dataset]][[method]] = do.call(rbind,tmp_lst) %>% as.data.frame %>% rename(V1 = "rank")
    ranking_LRgenes_lst_plot[[dataset]][[method]] %<>% mutate(FC_nSenderCells = params_grid$FC_nSenderCells, 
                                                              FC_nReceiverCells = params_grid$FC_nReceiverCells,
                                                                                 l = params_grid$l,
                                                                                 method = method,
                                                                                 dataset = dataset)
  }

  #### Combine all methods together
  # for precision recall plots
  statistics_results_lst_recallprecision_plot[[dataset]] = do.call(rbind, statistics_results_lst_recallprecision_plot[[dataset]])
  
  # for ranking plots
  ranking_LRgenes_lst_plot[[dataset]] = do.call(rbind, ranking_LRgenes_lst_plot[[dataset]])
  
  
  ###########################
  ##### Plot showing FP #####
  
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  FP = ggplot(tmp_df, aes(y = FP, x = method, color = method)) +
    geom_point() +
    scale_y_log10() +
    ylab("FP")
  
  
  #######################
  ##### UpSet plots #####
  
  # TODO
  master_lst_precision_recall$significant_interactions$cellchat
  
  ##########################################
  ##### Generate Precision recall plot #####
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]] %>%
    mutate(ratio_ReceiverSender = (FC_nReceiverCells/FC_nSenderCells) %>% log2) %>% 
    group_by(ratio_ReceiverSender,method ) %>%
    summarise_at(vars(precision,recall), mean)
  
  # dont average across ratio_ReceiverSender
  tmp_df2 = statistics_results_lst_recallprecision_plot[[dataset]] %>%
    mutate(ratio_ReceiverSender = (FC_nReceiverCells/FC_nSenderCells) %>% log2)
  
  p = ggplot(tmp_df2, aes(x = precision, y = recall)) +
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
    group_by(method,l ) %>%
    summarise_at(vars(f1score), mean) %>%
    mutate(f1score = round(f1score,3))
  
  heatmap1 = ggplot(d, aes(method, l, fill= f1score)) +
    geom_tile(color = "black") +
    geom_text(aes(label = f1score), color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("heatmap of f1score according to index l (the value of l param can be seen in config.yaml)")+
    scale_y_continuous(breaks=l_param_index) +
    ylab("Index of l parameter")
  
  p1 = ggplot(tmp_df, aes(x = pct_signal_added, y = f1score, color = method)) + 
    geom_point(size = 0.5) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ l)+
    ggtitle("F1score as a function of percentage of cells expressing L/R faceted by L param") 
  
  ##########################################
  ##### Generate ranking LR genes plot #####
  
  tmp_df = ranking_LRgenes_lst_plot[[dataset]] %>% 
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
    group_by(method,l) %>%
    summarise_at(vars(rank), list(rank = mean)) %>%
    mutate(rank = round(rank,2))
  
  # Heatmap of rank according to l radius
  heatmap2 = ggplot(tmp_df2, aes(method, l, fill= rank)) +
    geom_tile(color = "black") +
    geom_text(aes(label = rank), color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("rank of how each method retrieves inflated LR")+
    scale_y_continuous(breaks=l_param_index) +
    ylab("Index of l parameter") + 
    labs(fill = "log10(rank)")
  
  
  ##### Save plots
  pdf(file.path(path_results_dir ,paste0(dataset, "_recall_precision_plots.pdf")), width = 12, height = 7)
  p %>% print
  p1 %>% print
  heatmap1 %>% print
  heatmap2 %>% print
  p2 %>% print
  FP %>% print
  
  ###################################################################
  ##### summary spider charts across all parameter combinations #####
  ###################################################################
  df = statistics_results_lst_recallprecision_plot[[dataset]] %>%
    group_by(method) %>%
    summarise_at(vars(precision,recall,f1score), mean) %>%
    cbind(.,ranking_LRgenes_lst_plot[[dataset]] %>%
            group_by(method) %>%
            mutate(log10_rank = round(log10(rank),2)) %>%
            summarise_at(vars(log10_rank), mean)) %>%
    tibble::column_to_rownames(., var = "method") %>%
    select(-method)
  
  # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
  df = rbind(c(0,0,0,0), df) # minimum values
  df = rbind(c(1,1,1,5), df) # maximum values
  
  # Solid line colors (no transparency)
  line_colors = c(
    rgb(228, 26, 28, maxColorValue = 255),  # red
    rgb(55, 126, 184, maxColorValue = 255), # blue
    rgb(77, 175, 74, maxColorValue = 255),  # green
    rgb(152, 78, 163, maxColorValue = 255), # purple
    rgb(255, 127, 0, maxColorValue = 255),  # orange
    rgb(255, 255, 51, maxColorValue = 255), # yellow
    rgb(166, 86, 40, maxColorValue = 255)   # brown
  )
  
  # Transparent fill colors (alpha = 80 out of 255)
  fill_colors = c(
    rgb(228, 26, 28, alpha = 80, maxColorValue = 255),  # red
    rgb(55, 126, 184, alpha = 80, maxColorValue = 255), # blue
    rgb(77, 175, 74, alpha = 80, maxColorValue = 255),  # green
    rgb(152, 78, 163, alpha = 80, maxColorValue = 255), # purple
    rgb(255, 127, 0, alpha = 80, maxColorValue = 255),  # orange
    rgb(255, 255, 51, alpha = 80, maxColorValue = 255), # yellow
    rgb(166, 86, 40, alpha = 80, maxColorValue = 255)   # brown
  )
  
  
  radarchart(df,  axistype=1,
             pcol=line_colors , pfcol=fill_colors , plwd=2 , plty=1,
             vlcex=1.5 ,cglty=2,cglcol = "#726666",
             title = paste("Averaged metrics across all combination of parameters for dataset -",dataset ),
             caxislabels = rep("", 5))
  # Add a legend
  legend(x=1.5, y=0.75, legend = rownames(df[c(-1,-2),]), bty = "n", pch = 20,text.col = "black", 
         col=fill_colors,cex=1.25, pt.cex=4, text.width = 0.1)
  
  text(x = c(0.2,0.4,0.6,0.8,1), y = c(-0.05), labels = c(1,2,3,4,5), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  text(x = c(-0.2,-0.4,-0.6,-0.8,-1), y = c(-0.05), labels = c(0,0.25,0.5,0.75,1), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  text(x = c(-0.1), y = c(-0.2,-0.4,-0.6,-0.8,-1), labels = c(0,0.25,0.5,0.75,1), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  text(x = c(-0.1), y = c(0.2,0.4,0.6,0.8,1), labels = c(0,0.25,0.5,0.75,1), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  
  ###################################################################
  ##### summary spider charts across all parameter combinations #####
  ###################################################################
  df = statistics_results_lst_recallprecision_plot[[dataset]] %>%
    group_by(method) %>%
    summarise_at(vars(precision,recall,f1score), mean) %>%
    cbind(.,ranking_LRgenes_lst_plot[[dataset]] %>%
            group_by(method) %>%
            mutate(log10_rank = round(log10(rank),2)) %>%
            summarise_at(vars(log10_rank), mean)) %>%
    tibble::column_to_rownames(., var = "method") %>%
    select(-method) %>% 
    t
  
  # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
  df = rbind(c(0,0,0,0,0,0,0), df) # minimum values
  df = rbind(c(1,1,1,1,1,1,1), df) # maximum values
  
  # Solid line colors (no transparency)
  line_colors = c(
    rgb(228, 26, 28, maxColorValue = 255)  # red
  )
  
  # Transparent fill colors (alpha = 80 out of 255)
  fill_colors = c(
    rgb(228, 26, 28, alpha = 80, maxColorValue = 255)  # red
  )
  
  # PRECISION
  radarchart(df[c(1,2,3),] %>% as.data.frame,  axistype=1,
             pcol=line_colors , pfcol=fill_colors , plwd=2 , plty=1,
             vlcex=1.5 ,cglty=2,cglcol = "#726666",
             title = paste("Averaged PRECISION across all combination of parameters"),
             caxislabels = rep("", 5))
  
  text(x = c(-0.05), y = c(-0.2,-0.4,-0.6,-0.8,-1), labels = c(0,0.25,0.5,0.75,1), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  text(x = c(-0.05), y = c(0.2,0.4,0.6,0.8,1), labels = c(0,0.25,0.5,0.75,1), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  
  # RECALL
  radarchart(df[c(1,2,4),] %>% as.data.frame,  axistype=1,
             pcol=line_colors , pfcol=fill_colors , plwd=2 , plty=1,
             vlcex=1.5 ,cglty=2,cglcol = "#726666",
             title = paste("Averaged RECALL across all combination of parameters"),
             caxislabels = rep("", 5))
  
  text(x = c(-0.05), y = c(-0.2,-0.4,-0.6,-0.8,-1), labels = c(0,0.25,0.5,0.75,1), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  text(x = c(-0.05), y = c(0.2,0.4,0.6,0.8,1), labels = c(0,0.25,0.5,0.75,1), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  
  # RANK
  # Change min and max values for rank
  df[1,] = rep(5)
  df[2,] = rep(1)
  
  radarchart(df[c(1,2,6),] %>% as.data.frame,  axistype=1,
             pcol=line_colors , pfcol=fill_colors , plwd=2 , plty=1,
             vlcex=1.5 ,cglty=2,cglcol = "#726666",
             title = paste("Averaged log10(RANK) across all combination of parameters"),
             caxislabels = rep("", 5))
  
  text(x = c(-0.05), y = c(-0.2,-0.4,-0.6,-0.8,-1), labels = c(1,2,3,4,5), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  text(x = c(-0.05), y = c(0.2,0.4,0.6,0.8,1), labels = c(1,2,3,4,5), col = rgb(0, 0, 1, alpha = 0.5), cex = 1)
  dev.off()
  
}
