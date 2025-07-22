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
library(edgeR)
library(jsonlite)
library(SpatialExperiment)
library(ComplexHeatmap)
library(ggspavis)
source("scripts/helper_functions.R")

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
n_neighbors = config$semiSimulation$n_neighbors %>% unlist %>% as.integer

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
  params_grid = expand.grid(indexLR = indexLR, n_neighbors = n_neighbors)
  
  # iterate over the grid of parameters
  for(i in 1:nrow(params_grid))
  {
    x = params_grid[i,] %>% as.numeric ;  names(x) = c("indexLR","n_neighbors")
    
    naming = paste0("indexLR_",x["indexLR"],"_n_neighbors_",x["n_neighbors"])
    
    # load correct count file depending on the params_grid
    inflated_counts_file = inflated_counts_files %>% extract(grepl(paste0("n_neigbors_",x["n_neighbors"],"_indexLR_",x["indexLR"] , ".tsv$"), inflated_counts_files))
    inflated_counts = read.table(file.path(file_path,inflated_counts_file)) 
    
    # load correct simulated_cellmetadata file depending on the params_grid
    simulated_cellmetadata_file = simulated_cellmetadata_files %>% extract(grepl(paste0("n_neigbors_",x["n_neighbors"],"_indexLR_",x["indexLR"] , ".json$"), simulated_cellmetadata_files))
    master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]] = read_json(file.path(file_path,simulated_cellmetadata_file)) %>% convert_json_to_df
    
    # load correct simulated interactions file depending on the params_grid
    simulated_interactions_file = simulated_interactions_files %>% extract(grepl(paste0("n_neigbors_",x["n_neighbors"],"_indexLR_",x["indexLR"] , ".RDS$"), simulated_interactions_files))
    master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] = readRDS(file.path(file_path,simulated_interactions_file))
    
    master_lst_diagnosticPlots[[naming]][["CT1"]] = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]$CT1_CT2$ligand
    master_lst_diagnosticPlots[[naming]][["CT2"]] = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]$CT1_CT2$receptor
    
    ####################################################################
    ##### Generate density plots of simulated genes b/a simulation #####
    
    for(CT in c("CT1","CT2"))
    {
      # check if cellnames are ordered
      stopifnot(master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$cell_ID == colnames(inflated_counts))
      
      set.seed(1)
      gene_names = master_lst_diagnosticPlots[[naming]][[CT]] %>% unlist
      
      # only works if gene_names has no subunits. R1_R2 will not work correctly here
      for(gene in gene_names)
      {
        tmp_metadata = master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$metadata
        # create df with counts of inflated and original count matrices
        tmp_df = data.frame(inflated_counts = inflated_counts[gene,tmp_metadata$Celltype == CT] %>% as.numeric,
                            original_counts = original_counts[gene,tmp_metadata$Celltype == CT] %>% as.numeric) %>% melt
        
        # calculate mean value for both counts
        mu = ddply(tmp_df, "variable", summarise, grp.mean=mean(value))
        
        p = ggplot(tmp_df, aes(x=value, color=variable)) +
          geom_density()+
          geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),
                     linetype="dashed") +
          ggtitle(paste0("density of counts | lines - mean | n_neighbors:",x["n_neighbors"], "  gene:", gene)) +
          theme(plot.title = element_text(size = 7, face = "bold"),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.title=element_blank(),
                axis.text=element_text(size=6)) 
        
        diagnostic_gene_densityPlots[[dataset]][[naming]][[CT]][[gene]] = p
      }
      rm(tmp_metadata)
    }
    
    ###########################################################
    ##### Generate Plots of real FC after semi-simulation #####
    
    # Select file realFC after simulation
    realFC_aftersimulation_file = file_path %>% list.files(., pattern = "realFC_aftersimulation") %>% extract(grepl(paste0("n_neigbors_",x["n_neighbors"],"_indexLR_",x["indexLR"] , ".RDS$"), .))
    
    tmp_lst = readRDS(file.path(file_path,realFC_aftersimulation_file)) %>% 
      unlist %>% subset(.,!is.infinite(.)) %>% 
      plot_FCafter_semisimulation(. , indexLR = x["indexLR"] ,theoreticalFC = FC, n_neighbors = x["n_neighbors"])
    
    diagnostic_plots_realFC[[dataset]][[paste0("indexLR_",x["indexLR"])]][[paste0("n_neighbors_",x["n_neighbors"])]] = tmp_lst %>% pluck("plot")
    
    diagnostic_df_realFC[[dataset]][[paste0("index_",i)]] = tmp_lst[c(2,3,4)] %>% as.data.frame
    
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
  
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
                                                                      n_neighbors_param = n_neighbors, dataset = dataset, CT_toPlot = c("CT1","CT2"))
}


######################
##### Save plots #####
for(dataset in datasets)
{
  pdf(file.path(path_results_dir ,paste0(dataset, "_diagnostic_plots.pdf")), width = 12, height = 7)
  do.call(ggarrange,c(diagnostic_plots_perCT[[dataset]]$avelogcpm_fixed_n_neighbors, common.legend = TRUE)) %>% print
  
  x = 1:length(diagnostic_plots_realFC[[dataset]])
  sapply(x, function(x) {do.call(ggarrange,diagnostic_plots_realFC[[dataset]][[x]]) %>% print}) %>% print
  #diagnostic_plots_MeanVar[[dataset]] %>% print
  
  
  do.call(ggarrange,c(diagnostic_gene_densityPlots[[dataset]][[1]][["CT1"]],  
                      diagnostic_gene_densityPlots[[dataset]][[1]][["CT2"]], common.legend = TRUE)) %>% print
  
  #do.call(ggarrange,c(diagnostic_gene_densityPlots[[dataset]]$indexLR_2_n_neighbors_2[[1]],  
  #                    diagnostic_gene_densityPlots[[dataset]]$indexLR_2_n_neighbors_2[[2]], common.legend = TRUE)) %>% print
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
    metrics_files = file.path(path_output_dir,dataset,"metrics/",method) %>%
      list.files(., pattern = "f1_score")
    ranking_LRgenes_files = file.path(path_output_dir,dataset,"metrics/",method) %>% 
      list.files(., pattern = "ranking_LRgenes")
    
    master_lst_precision_recall = list()
    for(file in metrics_files)
    {
      master_lst_precision_recall[["metrics"]][[file]] = read.csv((file.path(path_output_dir,dataset,"metrics/",method,file))) %>% unlist
    }
    
    for(file in ranking_LRgenes_files)
    {
      master_lst_precision_recall[["ranking"]][[file]] = read.table((file.path(path_output_dir,dataset,"metrics/",method,file))) %>% unlist %>% as.numeric
    }
    
    
    
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] = do.call(rbind,master_lst_precision_recall[["metrics"]]) %>% as.data.frame
    ranking_LRgenes_lst_plot[[dataset]][[method]] = do.call(rbind,master_lst_precision_recall[["ranking"]]) %>% as.data.frame
    
    ################################################################################################
    ##### Average the metrics for different indexLR but across same combination of parameters  #####
    
    # extract rows that have the same parameters besides indexLR and then average the columns
    # generate the grid of parameters used for naming the list to generate output
    # this params_grid is different from previous because before we didnt look at methods and now we have different l_parameters for each method
    params_grid = expand.grid(n_neighbors = n_neighbors, l = l_param_index)
    
    tmp_lst = list()
    for(row in 1:nrow(params_grid))
    {
      x = params_grid[row,] %>% as.numeric ;  names(x) = c("n_neighbors","l")
      string = paste0("n_neigbors_",x["n_neighbors"],"_l_",x["l"])
      
      # extract filenames -> compute mean across columns and save to temporary list
      filenames = statistics_results_lst_recallprecision_plot[[dataset]][[method]] %>% 
        rownames() %>% 
        extract(grepl(string, .))
      tmp_lst[[row]] = statistics_results_lst_recallprecision_plot[[dataset]][[method]][filenames,] %>% colMeans()
      
      names(tmp_lst)[row] = string
    }
    
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] = do.call(rbind,tmp_lst) %>% as.data.frame
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] %<>% mutate(n_neighbors = params_grid$n_neighbors, 
                                                                                l = params_grid$l,
                                                                                method = method,
                                                                                dataset = dataset)
    ###############################
    ##### Ranking of LR genes #####
    
    tmp_lst = list()
    
    for(row in 1:nrow(params_grid))
    {
      x = params_grid[row,] %>% as.numeric ;  names(x) = c("n_neighbors","l")
      string = paste0("n_neigbors_",x["n_neighbors"],"_l_",x["l"])
      
      # extract filenames -> compute mean across columns and save to temporary list
      filenames = ranking_LRgenes_lst_plot[[dataset]][[method]] %>% 
        rownames() %>% 
        extract(grepl(string, .))
      tmp_lst[[row]] = ranking_LRgenes_lst_plot[[dataset]][[method]][filenames,] %>% mean
      
      names(tmp_lst)[row] = string
    }
    
    ranking_LRgenes_lst_plot[[dataset]][[method]] = do.call(rbind,tmp_lst) %>% as.data.frame %>% rename(V1 = "rank")
    ranking_LRgenes_lst_plot[[dataset]][[method]] %<>% mutate(n_neighbors = params_grid$n_neighbors, 
                                                                                 l = params_grid$l,
                                                                                 method = method,
                                                                                 dataset = dataset)
  }

  #### Combine all datasets together
  # for precision recall plots
  statistics_results_lst_recallprecision_plot[[dataset]] = do.call(rbind, statistics_results_lst_recallprecision_plot[[dataset]])
  
  # for ranking plots
  ranking_LRgenes_lst_plot[[dataset]] = do.call(rbind, ranking_LRgenes_lst_plot[[dataset]])
  
  ##########################################
  ##### Generate Precision recall plot #####
  
  # PLot precision recall curves when generating 1 plot per n_neighbors
  lst_precision_recall_by_n_neighbors_plots = list()
  for(n in unique(n_neighbors))
  {
    tmp_df = filter(statistics_results_lst_recallprecision_plot[[dataset]], n_neighbors == n)
    p = ggplot(tmp_df, aes(y = precision, x = recall , color = method)) + 
      geom_point(size = 1.5) + 
      geom_line() + 
      ggtitle(paste0("n_neighbors: ", n))+
      #geom_text_repel(aes(label = FC), size = 2.25,segment.linetype = 5,nudge_x = 0.005/max(tmp_df$recall)) + ggtitle(paste0("n_neighbors = ",n)) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))  + 
      theme(legend.title = element_text(size = 7), 
            legend.text = element_text(size = 7),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),  
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            plot.title = element_text(size=10))
    lst_precision_recall_by_n_neighbors_plots[[paste0("n_neighbors_",n)]] = p
  }
  
  ###################################
  ##### Generate f1 score plots #####
  
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  
  heatmap1 = ggplot(tmp_df, aes(method, l, fill= f1score)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("heatmap of f1score according to index l (the value of l param can be seen in config.yaml)")+
    scale_y_continuous(breaks=l_param_index) +
    ylab("Index of l parameter")
  
  heatmap2 = ggplot(tmp_df, aes(method, n_neighbors, fill= f1score)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("heatmap of f1score sorted by increasing values of n_neighbors")+
    scale_y_continuous(breaks=n_neighbors)
  
  p1 = ggplot(tmp_df, aes(x = n_neighbors, y = f1score, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~l) +
    ggtitle("Faceted by l")
  
  p2 = ggplot(tmp_df, aes(x = l, y = f1score, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~n_neighbors) +
    ggtitle("Faceted by n_neighbors")

  ##########################################
  ##### Generate ranking LR genes plot #####
  
  tmp_df = ranking_LRgenes_lst_plot[[dataset]] %>% mutate(rank = log10(rank)+1)
  
  p3 = ggplot(tmp_df %>% group_by(method,n_neighbors) %>% dplyr::summarise(rank = mean(rank))
              , aes(x = n_neighbors, y = rank, color = method)) + 
    geom_point() +
    geom_line() +
    ylab("log10(rank) + 1") +
    ggtitle("rank of inflated LR pair in the methods output")
  
  p4 = ggplot(tmp_df %>% group_by(method,l) %>% dplyr::summarise(rank = mean(rank))
              , aes(x = l, y = rank, color = method)) + 
    geom_point() +
    geom_line() +
    ylab("log10(rank) + 1")+
    ggtitle("rank of inflated LR pair in the methods output")
  
  heatmap3 = ggplot(tmp_df, aes(method, l, fill= rank)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("rank of how each method retrieves inflated LR")+
    scale_y_continuous(breaks=l_param_index) +
    ylab("Index of l parameter") + 
    labs(fill = "log10(rank) +1")
  
  heatmap4 = ggplot(tmp_df, aes(method, n_neighbors, fill= rank)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("rank of how each method retrieves inflated LR")+
    scale_y_continuous(breaks=l_param_index) +
    ylab("n_neighbors parameter")+ 
    labs(fill = "log10(rank) +1")
  
  ##### Save plots
  pdf(file.path(path_results_dir ,paste0(dataset, "_recall_precision_plots.pdf")), width = 12, height = 7)
  ggarrange(plotlist = lst_precision_recall_by_n_neighbors_plots, common.legend = T) %>% print
  p1 %>% print
  p2 %>% print
  p3 %>% print
  p4 %>% print
  heatmap1 %>% print
  heatmap2 %>% print
  heatmap3 %>% print
  heatmap4 %>% print
  dev.off()
}

