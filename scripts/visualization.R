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
FC = config$semiSimulation$FC %>% unlist %>% as.double
vec_n_neighbors = config$semiSimulati$n_neighbors %>% unlist %>% as.integer

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
  params_grid = expand.grid(FC = FC, n_neighbors = vec_n_neighbors)
  
  # iterate over the grid of parameters
  for(i in 1:nrow(params_grid))
  {
    x = params_grid[i,] %>% as.numeric ;  names(x) = c("FC","n_neighbors")
    
    naming = paste0("FC_",x["FC"],"_n_neighbors_",x["n_neighbors"])
    
    # load correct count file depending on the params_grid
    inflated_counts_file = inflated_counts_files %>% extract(grepl(paste0("_" , x["FC"] , "_", ".*",x["n_neighbors"] , ".tsv$"), inflated_counts_files))
    inflated_counts = read.table(file.path(file_path,inflated_counts_file)) 
    
    # load correct simulated_cellmetadata file depending on the params_grid
    simulated_cellmetadata_file = simulated_cellmetadata_files %>% extract(grepl(paste0("_" , x["FC"] , "_", ".*",x["n_neighbors"] , ".json$"), simulated_cellmetadata_files))
    master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]] = read_json(file.path(file_path,simulated_cellmetadata_file)) %>% convert_json_to_df
    
    # load correct simulated interactions file depending on the params_grid
    simulated_interactions_file = simulated_interactions_files %>% extract(grepl(paste0("_" , x["FC"] , "_",".*",x["n_neighbors"] , ".RDS$"), simulated_interactions_files))
    master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] = readRDS(file.path(file_path,simulated_interactions_file))
    
    ########################################
    ##### Generate L/R inflated per CT #####
    for(CT in c("CT1","CT2"))
    {
      # Find Ligand genes which were inflated in specified celltype 
      n = which(names(master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]) %>% str_split(.,"_") %>% lapply(.,"[[",1) %>% unlist %in% CT)
      L_genes_inflated_per_CT_to_keep = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]][n] %>% unlist %>% str_split(.,"_") %>% 
        lapply("[[", 1) %>% unlist %>% setdiff(., "subunit") # remove subunit string
      
      # Find Receptor genes which were inflated in specified celltype inflated 
      n = which(names(master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]) %>% str_split(.,"_") %>% lapply(.,"[[",2) %>% unlist %in% CT)
      R_genes_inflated_per_CT_to_keep = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]][n] %>% unlist %>% str_split(.,"_") %>% 
        lapply("[[", 2) %>% unlist%>% setdiff(., "subunit") # remove subunit string
      
      # Save inflated genes per CT
      master_lst_diagnosticPlots[[naming]][[CT]] = list(L = L_genes_inflated_per_CT_to_keep %>% unique , R = R_genes_inflated_per_CT_to_keep %>% unique)
    }
    
    ####################################################################
    ##### Generate density plots of simulated genes b/a simulation #####
    
    for(CT in c("CT1","CT2"))
    {
      # check if cellnames are ordered
      stopifnot(master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$cell_ID == colnames(inflated_counts))
      
      set.seed(1)
      gene_names = master_lst_diagnosticPlots[[naming]][[CT]] %>% unlist %>% sample(.,12) # sample only 9 genes as the report will contain only 9
      
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
          ggtitle(paste0("density of counts | lines - mean | FC:", x["FC"], "  n_neighbors:",x["n_neighbors"], "  gene:", gene)) +
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
    realFC_aftersimulation_file = file_path %>% list.files(., pattern = "realFC_aftersimulation") %>% extract(grepl(paste0("_",x["FC"] , "_" ,".*",x["n_neighbors"] , ".RDS$"), .))
    
    tmp_lst = readRDS(file.path(file_path,realFC_aftersimulation_file)) %>% 
      unlist %>% subset(.,!is.infinite(.)) %>% 
      plot_FCafter_semisimulation(. , theoreticalFC = x["FC"], n_neighbors = x["n_neighbors"])
    
    diagnostic_plots_realFC[[dataset]][[paste0("FC_",x["FC"])]][[paste0("n_neighbors_",x["n_neighbors"])]] = tmp_lst %>% pluck("plot")
    
    diagnostic_df_realFC[[dataset]][[paste0("index_",i)]] = tmp_lst[c(2,3,4)] %>% as.data.frame
    
    ##########################
    ##### Save AvelogCPM #####
    # filter the inflated counts based to contain cells belonging to the celltype indicated by CT_toPlot
    # This is for "compute_diagnostic_plots" function in order to save memory instead of saving entire inflated counts
    master_lst_diagnosticPlots[[naming]][["inflated_counts_aveLogCPM"]] = inflated_counts[,master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$metadata$Celltype %in% "CT1"] %>% aveLogCPM
    names(master_lst_diagnosticPlots[[naming]][["inflated_counts_aveLogCPM"]]) = rownames(inflated_counts)
  }
  
  #########################################################
  ##### Compute avelogCPM before and after simulation #####
  
  # filter original avelogcpm counts to contain same genes as the inflated count matrices
  original_counts = original_counts[which(rownames(original_counts) %in% rownames(inflated_counts)),]
  
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots, 
                                                               FC_param = FC, n_neighbors_param = vec_n_neighbors, dataset = dataset, 
                                                               metadata = master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$metadata, CT_toPlot = "CT1")
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
  
  
  do.call(ggarrange,c(diagnostic_gene_densityPlots[[dataset]]$FC_2_n_neighbors_6[[1]], common.legend = TRUE)) %>% print
  #do.call(ggarrange,c(diagnostic_gene_densityPlots[[dataset]]$FC_10_n_neighbors_4[[1]], common.legend = TRUE)) %>% print
  #do.call(ggarrange,c(diagnostic_gene_densityPlots[[dataset]]$FC_10_n_neighbors_4[[1]], common.legend = TRUE)) %>% print
  #do.call(ggarrange,c(diagnostic_gene_densityPlots[[dataset]]$FC_10_n_neighbors_4[[1]], common.legend = TRUE)) %>% print
  
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
    CT_statistics_files = file.path(path_output_dir,dataset,"metrics/",method) %>%
      list.files(., pattern = "CT_statistics")
    master_lst_precision_recall = list()
    for(file in CT_statistics_files)
    {
      master_lst_precision_recall[[file]] = read.csv((file.path(path_output_dir,dataset,"metrics/",method,file))) %>% unlist
    }
    
    #############################################################################
    ##### Ranking of LR genes across top n of significant hits from methods #####
    '
    tmp_ranking_LRgenes_lst = list()
    ranking_LRgenes_files = file.path(path_output_dir,dataset,"metrics/",method) %>% list.files(., pattern = "ranking_LRgenes")
    
    for(file in ranking_LRgenes_files)
    {
      tmp_ranking_LRgenes_lst[[file]] = read.csv((file.path(path_output_dir,dataset,"metrics/",method,file))) 
    }
    ranking_LRgenes_lst_plot[[dataset]][[method]] = do.call(rbind,tmp_ranking_LRgenes_lst) %>% mutate(. , method = method)
    
    rm(tmp_ranking_LRgenes_lst)
    '
    ##### replace the theoretical FC for the median of all effective FC for all genes
    tmp_diagnostic_df_realFC = do.call(rbind, diagnostic_df_realFC[[dataset]])
    
    # Add FC and n_neighbors columns
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] = do.call(rbind,master_lst_precision_recall) %>% as.data.frame
    
    x = rownames(statistics_results_lst_recallprecision_plot[[dataset]][[method]]) %>% stringr::str_split("_")
    tmp_FC = lapply(x, "[[" , 4) %>% unlist %>% as.double
    n_neighbors = lapply(x, "[[" , 7) %>% unlist %>% as.double
    
    # replace the theoretical FC for the median of effective FC across genes
    FC_real = vector()
    for(i in 1:length(n_neighbors))
    {
      FC_real = append(FC_real,filter(tmp_diagnostic_df_realFC, theoreticalFC == tmp_FC[i] & n_neighbors == n_neighbors[i]) %>% 
                         select("FC_real_median") %>% as.numeric) %>% round
    }
    rm(tmp_diagnostic_df_realFC)
    
    statistics_results_lst_recallprecision_plot[[dataset]][[method]] %<>% mutate(.,FC = FC_real, n_neighbors = n_neighbors, method = method)
  }

    # for precision recall plots
  statistics_results_lst_recallprecision_plot[[dataset]] = do.call(rbind, statistics_results_lst_recallprecision_plot[[dataset]]) %>% mutate(., dataset = dataset)
  
  # Add FC and n_neighbors info for ranking LR genes plot 
  '
  ranking_LRgenes_lst_plot[[dataset]] = ranking_LRgenes_lst_plot[[dataset]] %>% do.call(rbind,.) %>%
    mutate(n_neighbors = statistics_results_lst_recallprecision_plot[[dataset]]$n_neighbors ,
           FC = statistics_results_lst_recallprecision_plot[[dataset]]$FC)
  '
  ##########################################
  ##### Generate Precision recall plot #####
  
  n_neighbors = as.factor(statistics_results_lst_recallprecision_plot[[dataset]]$n_neighbors)
  FC = as.factor(statistics_results_lst_recallprecision_plot[[dataset]]$FC)
  
  # PLot precision recall curves when generating 1 plot per n_neighbors
  lst_precision_recall_by_n_neighbors = list()
  for(n in levels(n_neighbors))
  {
    data = filter(statistics_results_lst_recallprecision_plot[[dataset]], n_neighbors == n)
    p = ggplot(data, aes(y = precision, x = recall , color = method)) + 
      geom_point(size = 1.5) + 
      geom_line() + 
      geom_text_repel(aes(label = FC), size = 2.25,segment.linetype = 5,nudge_x = 0.005/max(data$recall)) + ggtitle(paste0("n_neighbors = ",n)) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01))  + 
      theme(legend.title = element_text(size = 7), 
            legend.text = element_text(size = 7),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),  
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            plot.title = element_text(size=10))
    lst_precision_recall_by_n_neighbors[[paste0("n_neighbors_",n)]] = p
  }
  
  ##################################
  ##### Generate f1 score plot #####
  
  tmp_df = statistics_results_lst_recallprecision_plot[[dataset]]
  p1 = ggplot(tmp_df, aes(x = n_neighbors, y = f1score, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~FC + n_neighbors) +
    ggtitle("Faceted by FC")
  
  p2 = ggplot(tmp_df, aes(x = FC, y = f1score, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~n_neighbors) +
    ggtitle("Faceted by n_neighbors")
  
  ##########################################
  ##### Generate ranking LR genes plot #####
  '
  tmp_df = ranking_LRgenes_lst_plot[[dataset]]
  p3 = ggplot(tmp_df, aes(x = FC, y = ratio, color = method)) + 
    geom_point() +
    geom_line() +
    facet_grid(~n_neighbors) +
    ggtitle("Faceted by n_neighbors . Ratio - n_simulated_LR in top50_significant_LR")
  '
  ##### Save plots
  pdf(file.path(path_results_dir ,paste0(dataset, "_recall_precision_plots.pdf")), width = 12, height = 7)
  ggarrange(plotlist = lst_precision_recall_by_n_neighbors, common.legend = T) %>% print
  p1 %>% print
  p2 %>% print
  #p3 %>% print
  dev.off()
}
