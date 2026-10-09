# ============================================================================
# Final step of the pipeline. Reads output/final_scores.RDS (produced by
# metric_f1score_rankingLRgenes.R) plus the intermediate semi-simulation and
# method output files, and produces all diagnostic/summary figures:
#   - Diagnostic plots: AveLogCPM and per-gene density before/after semi-
#     simulation, to sanity-check how much signal was actually added.
#   - Precision/recall/F1 plots, faceted by method, l parameter, and
#     PCE_Sender/PCE_Receiver.
#   - UpSet plots of which LR interactions each method calls significant.
#   - Heatmaps of F1/recall across PCE_Sender x PCE_Receiver x l parameter,
#     and of how many Receiver cells / how much of the neighborhood each
#     method "sees".
#   - A cross-dataset summary plot, and a comparison of the spatial scattering vs. Colocalization semi-simulation strategies.
# ============================================================================

suppressMessages({
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

# ==============================================================================
# LOCAL HELPER FUNCTIONS
# ==============================================================================

# Used for upset1 / upset2
make_upset_plot = function(binary_matrix, column_title = "All significant interactions retrieved across all combination of parameters")
{
  m = make_comb_mat(binary_matrix)
  UpSet(
    m,
    top_annotation = upset_top_annotation(
      m,
      height = unit(12.5, "cm"),
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 12),
      gp = gpar(fill = "steelblue")
    ),
    pt_size = unit(3, "mm"),
    lwd = 1,
    right_annotation = upset_right_annotation(
      m,
      width = unit(3, "cm"),
      gp = gpar(fill = "darkred")
    ),
    row_names_gp = gpar(fontsize = 12),
    row_gap = unit(0, "mm"),
    column_title = column_title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold")
  )
}

# Used for heatmap2 / heatmap3 / heatmap4
make_pce_heatmap = function(matrix_data, name, col_split, top_ann, target_width_mm, target_height_mm, adaptive_lwd = FALSE, ...)
{
  n_cols = ncol(matrix_data)
  n_rows = nrow(matrix_data)
  dynamic_cell_width  = target_width_mm / n_cols
  dynamic_cell_height = target_height_mm / n_rows
  
  rect_lwd = if (adaptive_lwd) { if (dynamic_cell_height < 3) 0 else 1 } else 1
  
  ht = Heatmap(matrix_data,
               name = name,
               top_annotation = top_ann,
               column_split = col_split,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               width = n_cols * unit(dynamic_cell_width, "mm"),
               height = n_rows * unit(dynamic_cell_height, "mm"),
               row_names_side = "left",
               column_names_gp = gpar(fontsize = 10),
               show_column_names = FALSE,
               rect_gp = gpar(col = "white", lwd = rect_lwd))
  
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "bottom",
                       merge_legends = FALSE, ...)
}

# Used for p3 / p4 / p5 / p6 
make_f1_vs_coverage_plot = function(df, x_var, extra_styling = FALSE)
{
  p = ggplot(df, aes(x = .data[[x_var]], y = f1score, color = method, linetype = AcountsSpatialDistance)) +
    geom_point(size = 0.25) +
    geom_smooth(se = FALSE, span = 0.5) +
    facet_wrap(~ radius_param_index,
               labeller = as_labeller(scenario_labels)) +
    ggtitle("F1score faceted by radius") +
    scale_linetype_manual(values = c("FALSE" = "dashed", "TRUE" = "solid")) +
    theme_bw()
  
  if (extra_styling) {
    p = p +
      theme(legend.key.width = unit(1.5, "cm")) +
      guides(linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
      scale_color_manual(values = method_cols, name = "Methods")
  }
  return(p)
}

# Used for heatmap2_faceted_f1score / heatmap2_faceted_recall 
make_faceted_pce_heatmap = function(df, fill_var, fill_label)
{
  ggplot(df, aes(x = factor(PCE_Sender), y = factor(PCE_Receiver), fill = .data[[fill_var]])) +
    geom_tile() +
    facet_grid(radius_param_index ~ method, labeller = labeller(radius_param_index = scenario_labels, method = label_both))
    scale_fill_viridis_c(option = "magma", limits = c(0, 1)) + # 'magma' is great for seeing F1 hotspots
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey90"),
      panel.spacing = unit(0.5, "lines")
    ) +
    labs(
      title = "Comparing methods across radius and PCE_Sender/PCE_Receiver thresholds",
      x = "PCE_Sender",
      y = "PCE_Receiver",
      fill = fill_label
    )
}

# Used for indexLR_plot1 / indexLR_plot2
make_indexLR_ridge_plot = function(df, x_var, x_label)
{
  ggplot(df, aes(x = .data[[x_var]], y = filename_without_indexLR, fill = filename_without_indexLR)) +
    geom_density_ridges(stat = "binline", bins = 30, scale = 0.9, alpha = 0.85, draw_baseline = FALSE) +
    facet_wrap(~method, ncol = 2) +
    scale_x_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
    theme_bw() +
    labs(
      title = "Method performance across different indexLR",
      x = x_label,
      y = "PC Snd/Rcv parameters"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# ==============================================================================
# STEP 0: load config, params, and the metric results produced upstream
# ==============================================================================

option_list = list(
  make_option("--config.yaml", type = "character", default = "config.yaml", help = "Path to config.yaml"),
  make_option("--strategy", type = "character", help = "spatialScattering or spatialColocalization")
)
opt = parse_args(OptionParser(option_list = option_list))

strategy = opt$strategy

path_output_dir = file.path("output", strategy)
path_results_dir = file.path("output", strategy, "results")
dir.create(path_results_dir, recursive = TRUE, showWarnings = FALSE)
metric_results = readRDS(file.path(path_output_dir, "final_scores.RDS"))

config = yaml::read_yaml("config.yaml")

# set colors
method_cols = setNames(scales::hue_pal()(8),
                       c("cellchat", "cellphonedbv5", "lianaP_morans", "mistyR",
                         "NICHES", "seurat_wilcoxon", "spatialdm", "stlearn"))

# replace nomenclature for radius with scenarions
scenario_labels = c(`0` = "negative", `1` = "optimal", `2` = "extended")

methods = config$methods %>% unlist
datasets = config$datasets %>% unlist
indexLR = config$indexLR_toSample %>% unlist %>% as.double
radius_param_index = config$radius_param_index %>% unlist %>% as.double
PCE_Sender = config$semiSimulation$PCE_Sender %>% unlist %>% as.double
PCE_Receiver = config$semiSimulation$PCE_Receiver %>% unlist %>% as.double

# ==============================================================================
# STEP 1: diagnostic plots (AveLogCPM + per-gene density, before/after
# semi-simulation), one PDF per dataset
# ==============================================================================

master_lst_diagnosticPlots = list()
diagnostic_gene_densityPlots = list()
diagnostic_plots_real = list()
diagnostic_plots_perCT = list()
gene_metadata_lst = list()

for (dataset in datasets)
{
  # Select inflated count files
  file_path = file.path(path_output_dir, paste0(dataset, "_semiSimulation_NB"))
  inflated_counts_files = file_path %>% list.files(., pattern = "inflated_counts")
  
  # Select all simulated interactions files
  simulated_interactions_files = list.files(file_path, pattern = "simulated_interactions")
  
  # Select simulated cellmetadata files
  cellmetadata_path = file.path(path_output_dir, paste0(dataset, "_semiSimulation_NB"))
  simulated_cellmetadata_files = list.files(cellmetadata_path, pattern = "simulated_cellmetadata")
  
  # Read original (pre-simulation) counts
  original_counts = read.table(file.path("data/processed", strategy,dataset, paste0("processed_counts_", dataset, ".tsv")))
  
  # generate the grid of parameters used for naming the list to generate outputs
  params_grid = expand.grid(indexLR = indexLR, PCE_Sender = PCE_Sender, PCE_Receiver = PCE_Receiver)
  
  # Only build the per-gene diagnostic density plots (below) for a random
  # handful of parameter combinations
  n_density_plot_combos = min(6, nrow(params_grid))
  density_plot_combo_indices = sample(nrow(params_grid), n_density_plot_combos)
  
  for (i in 1:nrow(params_grid))
  {
    x = params_grid[i, ] %>% as.numeric
    names(x) = c("indexLR", "PCE_Sender", "PCE_Receiver")
    
    naming = paste0("PCE_Sender_", x["PCE_Sender"], "_PCE_Receiver_", x["PCE_Receiver"], "_indexLR_", x["indexLR"])
    
    # load the count/cellmetadata/simulated-interactions files matching this parameter combination
    inflated_counts_file = inflated_counts_files %>% magrittr::extract(grepl(paste0("_PCE_Sender_", x["PCE_Sender"], "_PCE_Receiver_", x["PCE_Receiver"], "_indexLR_", x["indexLR"], ".tsv$"), inflated_counts_files))
    inflated_counts = read.table(file.path(file_path, inflated_counts_file))
    
    simulated_cellmetadata_file = simulated_cellmetadata_files %>% magrittr::extract(grepl(paste0("_PCE_Sender_", x["PCE_Sender"], "_PCE_Receiver_", x["PCE_Receiver"], "_indexLR_", x["indexLR"], ".json$"), simulated_cellmetadata_files))
    master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]] = read_json(file.path(file_path, simulated_cellmetadata_file)) %>% convert_json_to_df
    
    simulated_interactions_file = simulated_interactions_files %>% magrittr::extract(grepl(paste0("_PCE_Sender_", x["PCE_Sender"], "_PCE_Receiver_", x["PCE_Receiver"], "_indexLR_", x["indexLR"], ".RDS$"), simulated_interactions_files))
    master_lst_diagnosticPlots[[naming]][["simulated_interactions"]] = readRDS(file.path(file_path, simulated_interactions_file))
    
    master_lst_diagnosticPlots[[naming]][["Sender"]] = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]$Sender_Receiver$ligand
    master_lst_diagnosticPlots[[naming]][["Receiver"]] = master_lst_diagnosticPlots[[naming]][["simulated_interactions"]]$Sender_Receiver$receptor
    
    # ------------------------------------------------------------------
    # generate density plots of simulated genes before/after simulation
    # ------------------------------------------------------------------
    for (CT in c("Sender", "Receiver"))
    {
      # check if cellnames are ordered
      stopifnot(master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$cell_ID == colnames(inflated_counts))
      
      set.seed(1)
      gene_names = master_lst_diagnosticPlots[[naming]][[CT]] %>% unlist
      
      # in case there are subunits, separate the genes
      if (grepl("_", gene_names)) {
        gene_names = str_split(gene_names, "_") %>% unlist
      }
      
      for (gene in gene_names)
      {
        # cumulative expression for cells where we added signal and all other cells
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
        
        tmp_metadata = master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$metadata
        
        if (i %in% density_plot_combo_indices)
        {
          # create df with counts of inflated and original count matrices
          # since we have Sender and Sender_signalAdded, use grepl to find both
          tmp_df = data.frame(inflated_counts = inflated_counts[gene, grepl(CT, tmp_metadata$Celltype)] %>% as.numeric,
                              original_counts = original_counts[gene, grepl(CT, tmp_metadata$Celltype)] %>% as.numeric) %>% melt
          
          mu = ddply(tmp_df, "variable", summarise, grp.mean = mean(value))
          
          p = ggplot(tmp_df, aes(x = value, color = variable)) +
            geom_density() +
            geom_vline(data = mu, aes(xintercept = grp.mean, color = variable),
                       linetype = "dashed") +
            ggtitle(gene, subtitle = paste0("PCE_S: ", x["PCE_Sender"], " | PCE_R: ", x["PCE_Receiver"])) +
            theme_light() +
            theme(
              plot.title   = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.title = element_blank(),
              legend.text  = element_text(size = 12),
              axis.text    = element_text(size = 10))
          
          diagnostic_gene_densityPlots[[dataset]][[naming]][[CT]][[gene]] = p
        }
      }
      rm(tmp_metadata)
    }
    
    # Save AveLogCPM (for compute_diagnostic_plots() below, instead of saving the whole inflated count matrix)
    master_lst_diagnosticPlots[[naming]][["inflated_counts_aveLogCPM"]] = inflated_counts[, master_lst_diagnosticPlots[[naming]][["simulated_cellmetadata"]]$metadata$Celltype %in% c("Sender", "Receiver")] %>% aveLogCPM
    names(master_lst_diagnosticPlots[[naming]][["inflated_counts_aveLogCPM"]]) = rownames(inflated_counts)
  }
  
  # Plot AveLogCPM before and after simulation
  # filter original avelogcpm counts to contain same genes as the inflated count matrices
  original_counts = original_counts[which(rownames(original_counts) %in% rownames(inflated_counts)), ]
  
  # generate plots for few combination of PCE_Sender/PCE_Receiver
  PCE_Sender_targets = sample(PCE_Sender, min(2, length(PCE_Sender)))
  PCE_Receiver_targets = sample(PCE_Receiver, min(2, length(PCE_Receiver)))
  
  # any one parameter combination's cellmetadata works here - base Sender/
  # Receiver Celltype assignment doesn't vary across combinations for a
  # given dataset/strategy
  representative_metadata = master_lst_diagnosticPlots[[1]]$simulated_cellmetadata$metadata
  
  diagnostic_plots_perCT[[dataset]] = compute_diagnostic_plots(counts = original_counts, master_lst = master_lst_diagnosticPlots,
                                                               metadata = representative_metadata,
                                                               PCE_Sender_targets = PCE_Sender_targets, PCE_Receiver_targets = PCE_Receiver_targets,
                                                               dataset = dataset, CT_toPlot = c("Sender", "Receiver"))
  
  
  # save plots
  pdf(file.path(path_results_dir, paste0(dataset, "_diagnostic_plots.pdf")), width = 12, height = 7)
  do.call(ggarrange, c(diagnostic_plots_perCT[[dataset]]$avelogcpm, common.legend = TRUE)) %>% print
  
  # generate density plots
  N = length(diagnostic_gene_densityPlots[[dataset]])
  sample_indices = sample(N, min(2, N))
  
  do.call(ggarrange, c(diagnostic_gene_densityPlots[[dataset]][[sample_indices[1]]][["Sender"]],
                       diagnostic_gene_densityPlots[[dataset]][[sample_indices[1]]][["Receiver"]],
                       diagnostic_gene_densityPlots[[dataset]][[sample_indices[2]]][["Sender"]],
                       diagnostic_gene_densityPlots[[dataset]][[sample_indices[2]]][["Receiver"]], common.legend = TRUE)) %>% print
  
  dev.off()
  
}

saveRDS(master_lst_diagnosticPlots , file.path(path_results_dir, "master_lst_diagnosticPlots.RDS"))
saveRDS(gene_metadata_lst, file.path(path_results_dir, "gene_metadata_lst.RDS"))

# ==============================================================================
# STEP 2: precision/recall/F1 plots, per dataset
# ==============================================================================

# per-parameter-combination results, still split out by indexLR (kept for the
# per-indexLR precision/recall plots later in this script)
statistics_results_lst_recallprecision_non_indexLRaveraged = metric_results %>% lapply(., function(dataset) {
  lapply(dataset, function(method) {
    tmp_df = method %>% do.call(rbind, .) %>%
      # Create the new column by removing the suffix
      # The regex "_indexLR_\\d+$" targets "_indexLR_" and all digits at the end
      mutate(filename = rownames(.))
  })
})

# same as above, but averaged across indexLR (rows with identical parameters
# besides indexLR get collapsed into one, taking the mean of numeric columns)
statistics_results_lst_recallprecision = metric_results %>% lapply(., function(dataset) {
  lapply(dataset, function(method) {
    tmp_df = method %>% do.call(rbind, .) %>%
      mutate(file_without_indexLR = str_remove(rownames(.), "_indexLR_\\d+$")) %>%
      group_by(file_without_indexLR) %>%
      dplyr::summarise(
        across(
          everything(),
          ~ if (is.numeric(.x)) mean(.x, na.rm = TRUE) else .x[1]
        ),
        .groups = "drop"
      )
  })
})

master_lst_precision_recall = list()
simulated_interactions_lst = list()

for (dataset in datasets)
{
  for (method in methods)
  {
    # Load files for upset plot
    significant_interactions_files = file.path(path_output_dir, dataset, method) %>%
      list.files(., pattern = "significant_interactions")
    simulated_interactions_files = file.path(path_output_dir, paste0(dataset, "_semiSimulation_NB")) %>%
      list.files(., pattern = "simulated_interactions")
    
    for (file in significant_interactions_files)
    {
      # retrieve all significant interactions per method
      master_lst_precision_recall[["significant_interactions"]][[method]][[file]] = read.table(file.path(path_output_dir, dataset, method, file), header = T, sep = "\t") %>%
        mutate(significant = as.logical(significant)) %>%
        filter(significant) %>%
        select(ligand_receptor) %>%
        unlist %>%
        unname
      
      # retrieve all significant among the simulated interactions per method
      indexLR = str_extract(file, "(?<=indexLR_)\\d+")
      
      n = simulated_interactions_files %>% grep(paste0("indexLR_", indexLR), .) %>%
        magrittr::extract(1) # the simulated interaction only changes with indexLR, so take the first
      
      # this stores same things multiple times
      # can be improved but for now it doesnt take much space nor time
      simulated_interactions_lst[[method]][[file]] = readRDS(file.path(path_output_dir, paste0(dataset, "_semiSimulation_NB"), simulated_interactions_files[n])) %>%
        unlist %>%
        str_flatten("_")
      
      master_lst_precision_recall[["intersect_significant_simulated"]][[method]][[file]] = ifelse(simulated_interactions_lst[[method]][[file]] %in% master_lst_precision_recall[["significant_interactions"]][[method]][[file]],
                                                                                                  simulated_interactions_lst[[method]][[file]], 0)
    }
  }
  
  data_across_methods = statistics_results_lst_recallprecision[[dataset]]
  
  # ------------------------------------------------------------------
  # % of cells from cellmetadata that we added signal to
  # ------------------------------------------------------------------
  # Here, radius_param_index parameter is not relevant (here we use it for matching with master list) nor indexLR as that doesnt affect amount of cells. Thus, all files in simulated_cellmetadata_file should contain same %
  # load correct simulated_cellmetadata file depending on the params_grid | since all datasets have same values for % of cells signal added, which simulated_cellmetadata to load is not importnat
  
  params_grid = expand.grid(PCE_Sender = PCE_Sender, PCE_Receiver = PCE_Receiver, radius_param_index = radius_param_index)
  tmp_lst = list()
  for (row in 1:nrow(params_grid))
  {
    x = params_grid[row, ] %>% as.numeric
    names(x) = c("PCE_Sender", "PCE_Receiver", "radius_param_index")
    
    naming = paste0("PCE_Sender_", x["PCE_Sender"], "_PCE_Receiver_", x["PCE_Receiver"], "_l_", x["radius_param_index"])
    
    simulated_cellmetadata_file = simulated_cellmetadata_files %>% magrittr::extract(grepl(paste0("_PCE_Sender_", x["PCE_Sender"], "_PCE_Receiver_", x["PCE_Receiver"]), simulated_cellmetadata_files))
    
    y = read_json(file.path(path_output_dir, paste0(dataset, "_semiSimulation_NB"), simulated_cellmetadata_file[1])) %>% convert_json_to_df # take only 1 file
    y = y$metadata
    # calculate percentage of "Sender_signal_added" and "Receiver_signal_added" cells
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
  reordered_list = tmp_lst[data_across_methods[[1]]$file_without_indexLR]
  
  # make sure the entire reordered vector has same names as the dataframe column
  stopifnot(!any(reordered_list %>% names %>% is.na))
  
  # Combine all methods together for precision/recall plots; add pct_signal_added variable
  data_across_methods = do.call(rbind, data_across_methods) %>%
    mutate(pct_SenderCells_expressing_L =
             rep(reordered_list %>% do.call(rbind.data.frame, .) %>% filter(Celltype == "Sender") %>% pull(pct_signal_added), length(methods)),
           N_SenderCells_expressing_L =
             rep(reordered_list %>% do.call(rbind.data.frame, .) %>% filter(Celltype == "Sender") %>% pull(amount_signal_added), length(methods)),
           pct_ReceiverCells_expressing_R =
             rep(reordered_list %>% do.call(rbind.data.frame, .) %>% filter(Celltype == "Receiver") %>% pull(pct_signal_added), length(methods)),
           N_ReceiverCells_expressing_R =
             rep(reordered_list %>% do.call(rbind.data.frame, .) %>% filter(Celltype == "Receiver") %>% pull(amount_signal_added), length(methods))
    )
  
  # ------------------------------------------------------------------
  # UpSet plots: all interactions retrieved, and simulated-among-significant
  # ------------------------------------------------------------------
  significant_interactions_retrieved = lapply(master_lst_precision_recall[["significant_interactions"]], function(method) {
    all_strings = unique(unlist(master_lst_precision_recall$significant_interactions))
    binary_df = as.data.frame(lapply(method, function(x) {
      as.numeric(all_strings %in% x)
    }))
    
    data.frame(method = ifelse(rowSums(binary_df) > 0, 1, 0), row.names = all_strings)
  }) %>% do.call(cbind, .)
  
  colnames(significant_interactions_retrieved) = names(master_lst_precision_recall$significant_interactions)
  
  intersect_significant_simulated_interactions = lapply(master_lst_precision_recall[["intersect_significant_simulated"]], function(method) {
    all_strings = unique(unlist(simulated_interactions_lst))
    binary_df = as.data.frame(lapply(method, function(x) {
      as.numeric(all_strings %in% x)
    }))
    
    data.frame(method = ifelse(rowSums(binary_df) > 0, 1, 0), row.names = all_strings)
  }) %>% do.call(cbind, .)
  
  # remove 0
  intersect_significant_simulated_interactions = intersect_significant_simulated_interactions[rownames(intersect_significant_simulated_interactions) != "0", ]
  
  colnames(intersect_significant_simulated_interactions) = names(master_lst_precision_recall$intersect_significant_simulated)
  
  upset1 = make_upset_plot(significant_interactions_retrieved)
  upset2 = make_upset_plot(intersect_significant_simulated_interactions)
  
  # ------------------------------------------------------------------
  # precision/recall density plot
  # ------------------------------------------------------------------
  tmp_df = data_across_methods %>%
    mutate(ratio_ReceiverSender = (PCE_Receiver / PCE_Sender) %>% log2)
  
  p1 = ggplot(tmp_df, aes(x = precision, y = recall)) +
    geom_density_2d(
      aes(fill = after_stat(level)),
      contour_var = "ndensity",
      h = c(0.05, 0.05) # fine tude KDE bandwidth as some methods have very low variance and density estimation doesnt work
    ) +
    facet_wrap(~ method, ncol = 4) +
    labs(title = "Independent 2D Density per method")
  
  # ------------------------------------------------------------------
  # F1 score plots
  # ------------------------------------------------------------------
  tmp_df = data_across_methods %>%
    mutate(AcountsSpatialDistance = ifelse(method %in% c("cellphonedbv5", "seurat_wilcoxon"), FALSE, TRUE)) # generate column showing if the method is spatial distance dependent or not
  
  d = tmp_df %>%
    group_by(method, radius_param_index) %>%
    summarise_at(vars(f1score), mean) %>%
    mutate(f1score = round(f1score, 3))
  
  heatmap1 = ggplot(d, aes(method, radius_param_index, fill = f1score)) +
    geom_tile(color = "black") +
    geom_text(aes(label = f1score), color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggtitle("heatmap of f1score according to index l (the value of l param can be seen in config.yaml)") +
    scale_y_continuous(breaks = radius_param_index) +
    ylab("Index of l parameter")
  
  p3 = make_f1_vs_coverage_plot(tmp_df, "pct_SenderCells_expressing_L", extra_styling = TRUE)
  p4 = make_f1_vs_coverage_plot(tmp_df, "N_SenderCells_expressing_L", extra_styling = TRUE)
  p5 = make_f1_vs_coverage_plot(tmp_df, "pct_ReceiverCells_expressing_R", extra_styling = FALSE)
  p6 = make_f1_vs_coverage_plot(tmp_df, "N_ReceiverCells_expressing_R", extra_styling = FALSE)
  
  # ------------------------------------------------------------------
  # heatmaps across PCE_Sender x PCE_Receiver x l parameter
  # ------------------------------------------------------------------
  tmp_df = data_across_methods
  
  heatmap2_faceted_f1score = make_faceted_pce_heatmap(tmp_df, "f1score", "F1-Score")
  heatmap2_faceted_recall = make_faceted_pce_heatmap(tmp_df, "recall", "Recall")
  
  # Create a 'Combination' label for columns
  tmp_df %<>%
    mutate(Param_Comb = paste0("PCE_S:", PCE_Sender, "\nPCE_R:", PCE_Receiver, "\nradius:", radius_param_index)) %>%
    select(method, Param_Comb, recall) %>%
    pivot_wider(names_from = Param_Comb, values_from = recall) %>%
    tibble::column_to_rownames("method") %>%
    as.matrix()
  
  # reorder columns numerically by PCE_S then PCE_R
  pce_s_num = as.numeric(str_extract(colnames(tmp_df), "(?<=PCE_S:)[0-9.]+"))
  pce_r_num = as.numeric(str_extract(colnames(tmp_df), "(?<=PCE_R:)[0-9.]+"))
  tmp_df = tmp_df[, order(pce_s_num, pce_r_num)]
  
  # Extract metadata for annotations (the labels at the top).
  # This keeps the PCE_Sender and PCE_Receiver values linked to the columns
  col_meta = data.frame(colnames(tmp_df)) %>%
    separate(1, into = c("PCE_S", "PCE_R", "radius"), sep = "\n") %>%
    mutate(PCE_S = str_remove(PCE_S, "^PCE_S:"),
           PCE_R = str_remove(PCE_R, "^PCE_R:"))
  
  col_meta$PCE_S  = factor(col_meta$PCE_S,  levels = str_sort(unique(col_meta$PCE_S),  numeric = TRUE))
  col_meta$PCE_R  = factor(col_meta$PCE_R,  levels = str_sort(unique(col_meta$PCE_R),  numeric = TRUE))
  col_meta$radius = factor(col_meta$radius, levels = str_sort(unique(col_meta$radius), numeric = TRUE))
  
  
  unique_pce_s = unique(col_meta$PCE_S)
  unique_pce_r = unique(col_meta$PCE_R)
  unique_radius = unique(col_meta$radius)
  
  pce_s_cols = setNames(brewer.pal(length(unique_pce_s), "Set1"), unique_pce_s)
  #pce_s_cols = pce_s_cols[1:2]
  pce_r_cols = setNames(brewer.pal(length(unique_pce_r), "Set2"), unique_pce_r)
  #pce_r_cols = pce_r_cols[1:2]
  radius_cols = setNames(brewer.pal(length(unique_radius), "Set2"), unique_radius)
  #radius_cols = radius_cols[1:2]
  
  top_ann = HeatmapAnnotation(
    PCE_S = col_meta$PCE_S,
    PCE_R = col_meta$PCE_R,
    col = list(
      PCE_S = pce_s_cols,
      PCE_R = pce_r_cols,
      radius = radius_cols
    )
  )
  
  # f1 score color gradient (0 to 1) - currently unused, kept for reference
  col_f1 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  
  heatmap2 = make_pce_heatmap(tmp_df, name = "Recall", col_split = col_meta$radius, top_ann = top_ann,
                              target_width_mm = 300, target_height_mm = 30,
                              legend_grouping = "adjusted", ht_gap = unit(10, "mm"))
  
  # --- Plot how many overall cells are seen by method (must come after heatmap2, reuses top_ann/col_meta) ---
  tmp_df = data_across_methods
  
  tmp_df %<>%
    mutate(average_cells_perSender_seen_byMethod = as.numeric(average_cells_perSender_seen_byMethod),
           average_cells_perSender_seen_byMethod = ifelse(method %in% c("seurat_wilcoxon", "cellphonedbv5"), NA, average_cells_perSender_seen_byMethod),
           Param_Comb = paste0("PCE_S:", PCE_Sender, "\nPCE_R:", PCE_Receiver, "\nradius:", radius_param_index)) %>%
    select(method, Param_Comb, average_cells_perSender_seen_byMethod) %>%
    pivot_wider(names_from = Param_Comb, values_from = average_cells_perSender_seen_byMethod) %>%
    tibble::column_to_rownames("method") %>%
    as.matrix()
  
  heatmap3 = make_pce_heatmap(tmp_df, name = "Mean_neighborhood_size", col_split = col_meta$radius, top_ann = top_ann,
                              target_width_mm = 300, target_height_mm = 30,
                              legend_grouping = "adjusted", ht_gap = unit(10, "mm"))
  
  # --- Plot how many Receiver cells are captured by method (must come after heatmap2, reuses top_ann/col_meta) ---
  tmp_df = data_across_methods
  
  tmp_df %<>%
    mutate(ratio_Receiver_seen_byMethod = ifelse(method %in% c("seurat_wilcoxon", "cellphonedbv5"), NA, ratio_Receiver_seen_byMethod),
           Param_Comb = paste0("PCE_S:", PCE_Sender, "\nPCE_R:", PCE_Receiver, "\nradius:", radius_param_index)) %>%
    select(method, Param_Comb, ratio_Receiver_seen_byMethod) %>%
    pivot_wider(names_from = Param_Comb, values_from = ratio_Receiver_seen_byMethod) %>%
    tibble::column_to_rownames("method") %>%
    as.matrix()
  
  heatmap4 = make_pce_heatmap(tmp_df, name = "receiver_coverage_pct", col_split = col_meta$radius, top_ann = top_ann,
                              target_width_mm = 300, target_height_mm = 30,
                              legend_grouping = "adjusted", ht_gap = unit(10, "mm"))
  
  
  # ------------------------------------------------------------------
  # ranking LR genes plot
  # ------------------------------------------------------------------
  tmp_df = data_across_methods %>%
    mutate(average_cells_perSender_seen_byMethod = as.numeric(average_cells_perSender_seen_byMethod) %>% as.integer) %>%
    mutate(is_spatial_aware = !method %in% c("seurat_wilcoxon", "cellphonedbv5"))
  
  
  # Reshape data so FP and f1score are in a single column
  tmp_df_long = tmp_df %>%
    pivot_longer(
      cols = c(FP, f1score),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(metric = case_when(
      metric == "FP" ~ "Number of Interactions Retrieved",
      metric == "f1score" ~ "F1-Score"
    ))
  
  N_interactions_f1score_plot = ggplot(tmp_df_long, aes(x = average_cells_perSender_seen_byMethod, y = value, color = method)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_smooth(se = FALSE, span = 0.5, size = 0.8) +
    scale_color_manual(values = method_cols) +
    facet_grid(metric ~ radius_param_index, scales = "free_y", labeller = labeller(radius_param_index = scenario_labels)) +
    scale_y_log10() +
    labs(
      x = "Mean_neighborhood_size",
      y = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right"
  )

  # ridge plot of rankings and detection rate
  raw_grid = do.call(rbind, statistics_results_lst_recallprecision_non_indexLRaveraged[[dataset]])
  
  detection_summary = raw_grid %>%
    mutate(normalized_rank = (rank - 1) / (N - 1)) %>%
    group_by(method, radius_param_index) %>%
    dplyr::summarise(
      detection_rate = mean(!is.na(rank)))
  
  detection_rate_labels = detection_summary %>%
    mutate(label = paste0(round(100 * detection_rate), "% detected"))
  
  rank_breaks = 10^(0:ceiling(log10(max(raw_grid$rank, na.rm = TRUE))))
  ridge_plot = ggplot(tmp_df, aes(x = rank, y = method, fill = method, linetype = is_spatial_aware)) +
    geom_density_ridges() +
    geom_text(data = detection_rate_labels, aes(x = Inf, y = method, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = -0.6, size = 4.2, fontface = "bold", color = "grey20") +
    scale_x_log10(breaks = rank_breaks, labels = scales::label_number(big.mark = "")) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_fill_manual(values = method_cols) +
    theme_bw() +
    ggtitle("rank and detection rate for every method in retrieving inflated LR") +
    facet_wrap(~ radius_param_index, labeller = as_labeller(scenario_labels))
  
  # ------------------------------------------------------------------
  # precision/recall plots across indexLR
  #
  # The goal here is to understand how values of recall and precision change
  # depending on which ligand/receptor pair (indexLR) one chooses..
  # ------------------------------------------------------------------
  
  # For the final plot, too many combinations of parameters are used, reduce them here
  n = length(PCE_Sender)
  PCE_Sender_touse = PCE_Sender[unique(c(seq(1, n, by = 3), n))]
  n = length(PCE_Receiver)
  PCE_Receiver_touse = PCE_Receiver[unique(c(seq(1, n, by = 3), n))]
  
  # summarise data where we merge information of how many cells per celltype and
  # combination of parameters are expressing ligands/receptors, and the
  # f1score/precision/recall scores for further plotting
  lst_data = list()
  for (indexLR in config$indexLR_toSample)
  {
    # cell number statistics for Sender (expressing ligands)
    ligands = lapply(gene_metadata_lst[[dataset]], function(x) {x[["Sender"]] %>% do.call(rbind.data.frame, .)
    }) %>%
      do.call(rbind.data.frame, .) %>%
      tibble::rownames_to_column("filename") %>%
      mutate(filename = gsub("(_indexLR_).*", "", filename)) %>%
      filter(celltype == "Sender") %>%
      dplyr::group_by(indexLR, filename) %>%
      dplyr::summarise(
        across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
        across(where(is.character), ~paste(.x, collapse = ", "))
      ) %>% mutate(filename = str_c(filename, "_indexLR_", indexLR))
    
    # cell number statistics for Receiver (expressing receptors)
    receptors = lapply(gene_metadata_lst[[dataset]], function(x) {x[["Receiver"]] %>% do.call(rbind.data.frame, .)
    }) %>%
      do.call(rbind.data.frame, .) %>%
      tibble::rownames_to_column("filename") %>%
      mutate(filename = gsub("(_indexLR_).*", "", filename)) %>%
      filter(celltype == "Receiver") %>%
      dplyr::group_by(indexLR, filename) %>%
      dplyr::summarise(
        across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
        across(where(is.character), ~paste(.x, collapse = ", "))
      ) %>% mutate(filename = str_c(filename, "_indexLR_", indexLR))
    
    # add f1score/precision/recall values for each combination of PCE parameters, per method
    for (m in names(statistics_results_lst_recallprecision_non_indexLRaveraged[[dataset]]))
    {
      tmp_df = statistics_results_lst_recallprecision_non_indexLRaveraged[[dataset]][[m]] %>%
        select(precision, recall, PCE_Sender, PCE_Receiver, filename) %>%
        filter(PCE_Sender %in% PCE_Sender_touse & PCE_Receiver %in% PCE_Receiver_touse) %>%
        mutate(filename = str_remove(filename, "_l_\\d+")) %>%
        dplyr::group_by(filename) %>%
        dplyr::summarise(
          across(where(is.numeric) | is.logical, ~mean(.x, na.rm = TRUE))
        )
      
      tmp_df2 = left_join(ligands, receptors, by = "filename")
      tmp_df3 = left_join(tmp_df, tmp_df2, by = "filename")
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
  final_df$filename_without_indexLR %<>% gsub("PCE_Sender_", "S", .) %>% gsub("PCE_Receiver_", "R", .)
  
  # simplify the filename column for easier plotting
  final_df$filename_without_indexLR = str_replace(final_df$filename_without_indexLR, "(?<=_R)(\\d+\\.?\\d*)", function(m) {
    sprintf("%.1f", as.numeric(m))
  }) %>% str_replace(., "(?<=S)(\\d+\\.?\\d*)", function(m) {
    sprintf("%.1f", as.numeric(m))
  }) %>% str_replace_all(., "\\.0", "")
  
  # generate levels
  level_order = unique(final_df$filename_without_indexLR)
  level_order = level_order[order(
    as.numeric(str_extract(level_order, "(?<=S)\\d+")),
    as.numeric(str_extract(level_order, "(?<=R)\\d+"))
  )]
  final_df$filename_without_indexLR = factor(final_df$filename_without_indexLR, levels = level_order)
  
  indexLR_plot1 = make_indexLR_ridge_plot(final_df, "precision", "Precision")
  indexLR_plot2 = make_indexLR_ridge_plot(final_df, "recall", "Recall")
  
  # ------------------------------------------------------------------
  # save the precision/recall plots for this dataset
  # ------------------------------------------------------------------
  pdf(file.path(path_results_dir, paste0(dataset, "_recall_precision_plots.pdf")), width = 15, height = 7)
  p1 %>% print
  ggarrange(plotlist = list(p3, p4), common.legend = T) %>% print
  heatmap1 %>% print
  heatmap2_faceted_f1score %>% print
  heatmap2_faceted_recall %>% print
  ridge_plot %>% print
  N_interactions_f1score_plot %>% print
  upset1 %>% print
  upset2 %>% print
  indexLR_plot1 %>% print
  indexLR_plot2 %>% print
  dev.off()
  
  pdf(file.path(path_results_dir, paste0(dataset, "_recall_precision_plots2.pdf")), width = 20, height = 7)
  heatmap2 %>% print
  heatmap3 %>% print
  heatmap4 %>% print
  dev.off()
}

# ==============================================================================
# STEP 3: cross-dataset summary plot
# ==============================================================================

tmp_df_long = lapply(names(statistics_results_lst_recallprecision), function(dataset) {
  do.call(rbind.data.frame, statistics_results_lst_recallprecision[[dataset]]) %>%
    mutate(dataset = dataset)
}) %>%
  do.call(rbind.data.frame, .) %>%
  mutate(
    PCE_Receiver = as.numeric(as.character(PCE_Receiver)),
    PCE_Sender = as.numeric(as.character(PCE_Sender))
  ) %>%
  arrange(PCE_Receiver, PCE_Sender) %>%
  mutate(params = paste0("PCE_Receiver_", PCE_Receiver, "_PCE_Sender_", PCE_Sender, "_l_", radius_param_index)) %>%
  mutate(params = factor(params, levels = unique(params)))

plot_data = tmp_df_long %>%
  mutate(
    # Step 1: Extract the numerical digits directly out of the parameter string
    rcv_val = str_split_i(params, "_", 3) %>% as.numeric(),
    snd_val = str_split_i(params, "_", 6) %>% as.numeric(),
    
    # Step 2: Calculate the product feature for the continuous X-axis
    Snd_Rcv_Product = rcv_val * snd_val,
    
    # Step 3: Format faceting variables as clean factors
    dataset = factor(dataset),
    radius_param_index = factor(paste0("radius_param:", radius_param_index)),
    
    # Optional: Build your Boolean metadata flag for spatial awareness
    AccountsSpatialDistance = !method %in% c("cellphonedbv5", "seurat_wilcoxon")
  )

# Define manual dashed patterns for the baseline (spatial-unaware) methods
all_methods = unique(plot_data$method)
line_types  = setNames(rep("solid", length(all_methods)), all_methods)
line_types["cellphonedbv5"]    = "dashed"
line_types["seurat_wilcoxon"] = "dashed"

product_f1score_plot = ggplot(plot_data, aes(x = Snd_Rcv_Product, y = f1score,
                                             color = method,
                                             linetype = method,
                                             linewidth = AccountsSpatialDistance,
                                             group = method)) +
  
  geom_smooth(method = "loess",
              span = 0.75,
              se = FALSE,
              alpha = 0.85) +
  
  facet_grid(dataset ~ radius_param_index, scales = "free_x") +
  
  scale_color_viridis_d(option = "D", name = "Methods") +
  scale_linetype_manual(values = line_types, name = "Methods") +
  
  # Keep the uniform linewidths as configured
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

ggsave(filename = file.path(path_results_dir, "summary_across_datasets.png"), plot = product_f1score_plot, width = 200, height = 150, units = "mm")

# ==============================================================================
# STEP 4: comparison of the Scattering vs. Colocalization semi-simulation
# strategies.
#
# This script runs once per strategy (see --strategy), but this particular
# plot needs *both* strategies' results at once - so it's guarded to only
# execute during the spatialScattering. Uses the real per-strategy output paths directly
# (output/spatialScattering/final_scores.RDS,
# output/spatialColocalization/final_scores.RDS)
# ==============================================================================
if (strategy == "spatialScattering")
{
  # make plot for precision/recall for colocalization and scattering
  out = lapply(datasets, function(dataset) {
    message(dataset)
    
    # SCATTERING
    metric_results_scattering = readRDS("output/spatialScattering/final_scores.RDS")
    statistics_results_lst_recallprecision_scattering = metric_results_scattering %>% lapply(., function(dataset) {
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
              ~ if (is.numeric(.x)) mean(.x, na.rm = TRUE) else .x[1]
            ),
            .groups = "drop"
          )
      })
    })
    
    # COLOCALIZATION
    metric_results_colocalization = readRDS("output/spatialColocalization/final_scores.RDS")
    statistics_results_lst_recallprecision_colocalization = metric_results_colocalization %>% lapply(., function(dataset) {
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
              ~ if (is.numeric(.x)) mean(.x, na.rm = TRUE) else .x[1]
            ),
            .groups = "drop"
          )
      })
    })
    
    data_across_methods_scattering = statistics_results_lst_recallprecision_scattering[[dataset]] %>% do.call(rbind, .)
    data_across_methods_colocalization = statistics_results_lst_recallprecision_colocalization[[dataset]] %>% do.call(rbind, .)
    # Assuming both dataframes have columns: "method", "recall", and "precision"
    data_across_methods_scattering$condition = "Scattering"
    data_across_methods_colocalization$condition = "Colocalization"
    data_across_methods_scattering$dataset = dataset
    data_across_methods_colocalization$dataset = dataset
    
    combined_data = bind_rows(data_across_methods_scattering, data_across_methods_colocalization)
    
    long_data = combined_data %>%
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
  pooled_data = long_data_plot %>%
    # Force Metric to be exactly what ggplot expects
    dplyr::group_by(method, condition, Metric, PCE_Sender, PCE_Receiver, radius_param_index) %>%
    dplyr::summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  # 2. Run the simplified plot
  p = ggplot(pooled_data, aes(x = method, y = Value, color = condition)) +
    geom_boxplot(
      position = position_dodge(width = 0.7),
      width = 0.6,
      alpha = 0.25, 
      outlier.shape = 21,
      outlier.size = 1.5,
      linewidth = 0.6
    ) +
    
    # This matches the capital "Metric" column from above
    facet_wrap(~Metric, scales = "free_y", nrow = 2) +
    
    scale_color_manual(values = c("#4292C6", "#EF3B2C")) +
    labs(
      title = "Benchmark performance profile: Scattering vs Colocalization",
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
  
  # saved outside any single strategy's results/ folder, since this plot
  # compares across both strategies rather than belonging to just one
  dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = "output/results/summary_across_datasets_scattering_vs_colocalization.png", plot = p, width = 200, height = 150, units = "mm")
}
