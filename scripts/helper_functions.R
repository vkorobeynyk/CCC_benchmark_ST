# ============================================================================
# Semi-simulation framework: inflates ligand expression in a chosen set of
# Sender cells and receptor expression in a chosen set of Receiver cells, by
# sampling from a negative binomial distribution fit to the gene's observed
# mean/dispersion. This creates a "ground truth" communication signal that
# CCC methods should be able to recover.
#
# Arguments:
#   counts                      - gene x cell count matrix to inflate
#   simulated_interactions_lst  - list holding the sampled Sender_Receiver
#                                  ligand/receptor pair(s) for this run
#   fraction_cells_expressingL  - target fraction of Sender cells that should
#                                  express the ligand after inflation
#   fraction_cells_expressingR  - target fraction of Receiver cells that
#                                  should express the receptor after inflation
#   genemetadata                - per-gene mean/dispersion estimates (edgeR)
#   metadata                    - cell metadata (must contain a Celltype column)
#   df_neighbors                - Sender -> Receiver neighbor lookup table
#
# Returns a list with:
#   counts_inflated - the count matrix with inflated L/R expression
#   cell_info       - metadata of which cells received inflated signal
# ============================================================================
semi_simulate = function(counts, simulated_interactions_lst, fraction_cells_expressingR, fraction_cells_expressingL,
                         genemetadata, metadata, df_neighbors)
{
  counts_inflated = counts
  cell_info       = list()
  means_perCT     = genemetadata$mean
  
  # the ligand/receptor gene(s) to inflate for this run (subunits are split, e.g. "R1_R2" -> "R1","R2")
  L_sample = simulated_interactions_lst$Sender_Receiver$ligand %>% str_split("_") %>% unlist
  R_sample = simulated_interactions_lst$Sender_Receiver$receptor %>% str_split("_") %>% unlist
  
  # NOTE on loop order: we reuse the variable `cells_toAdd_signal` across
  # iterations, so Receiver must always be processed *after* Sender - the
  # Receiver branch below reads the Sender cells chosen in the previous pass.
  for (tmp_CT in c("CTsender", "CTreceiver"))
  {
    # ------------------------------------------------------------------
    # Step 1: decide which genes to inflate and which cells receive signal
    # ------------------------------------------------------------------
    if (tmp_CT == "CTsender")
    {
      genes_to_sample = L_sample
      CT              = "Sender"
      
      n_sender_cells      = (length(colnames(df_neighbors)) * fraction_cells_expressingL) %>% floor
      cells_toAdd_signal  = sample(colnames(df_neighbors), n_sender_cells)
      
    } else if (tmp_CT == "CTreceiver")
    {
      genes_to_sample = R_sample
      CT              = "Receiver"
      
      # Sample Receiver cells according to the ratio fraction_cells_expressingR / fraction_cells_expressingL
      ratio_R_to_L = fraction_cells_expressingR / fraction_cells_expressingL
      
      if (ratio_R_to_L == 1)
      {
        # same number of Receivers as Senders -> just take the real neighbors of the chosen Senders
        cells_toAdd_signal = df_neighbors[, cells_toAdd_signal] %>% unlist %>% unname
        
      } else if (ratio_R_to_L > 1)
      {
        # more Receiver cells than Sender cells needed
        n_receiver_cells = (length(colnames(df_neighbors)) * fraction_cells_expressingR) %>% floor
        
        # first, take the real neighbors of the Sender cells that already carry signal
        neighbors_of_signal_senders = df_neighbors[, cells_toAdd_signal] %>%
          na.omit %>%
          t %>%
          as.data.frame %>%
          sample_n(n_sender_cells)
        
        # then top up with neighbors of Sender cells that were NOT chosen to carry signal
        neighbors_of_other_senders = df_neighbors[, which(!colnames(df_neighbors) %in% cells_toAdd_signal)] %>%
          na.omit %>%
          t %>%
          as.data.frame %>%
          sample_n(n_receiver_cells - n_sender_cells)
        
        cells_toAdd_signal = rbind(neighbors_of_signal_senders, neighbors_of_other_senders) %>% unlist %>% unname
        
      } else if (ratio_R_to_L < 1)
      {
        # fewer Receiver cells than Sender cells needed -> subsample the Senders' real neighbors
        n_receiver_cells = (length(colnames(df_neighbors)) * fraction_cells_expressingR) %>% floor
        
        cells_toAdd_signal = df_neighbors[, cells_toAdd_signal] %>%
          na.omit %>%
          t %>%
          as.data.frame %>%
          sample_n(n_receiver_cells) %>%
          unlist %>% unname
      }
    }
    
    # save which cells received inflated signal (used later for
    # diagnostics and for computing how much of the neighborhood each method "sees")
    cell_info[["Sender_Receiver"]][["cells_signalAdded"]][[CT]] = cells_toAdd_signal
    
    # ------------------------------------------------------------------
    # Step 2: inflate expression, gene by gene, for the chosen cells
    # ------------------------------------------------------------------
    for (gene_sample in genes_to_sample)
    {
      # zero out this gene's expression across the whole celltype (Sender or Receiver)
      # before adding the new, sampled expression values
      counts_inflated[gene_sample, metadata$Celltype %>% grepl(CT, .)] = 0
      
      # look up this gene's celltype-specific rate and dispersion (estimated via edgeR)
      gene_rate       = means_perCT[grep(paste0("^", gene_sample, "$"), means_perCT$gene_names), ] %>% .[CT] %>% as.numeric()
      gene_dispersion = genemetadata$disp %>% subset(gene == gene_sample) %>% select(edgeR_dispersion) %>% as.numeric
      
      # sample from a negative binomial size = 1/dispersion is the shape
      cell_offsets = genemetadata$edgeR_offsets$edgeR_offset[match(cells_toAdd_signal, genemetadata$edgeR_offsets$cell_ID)]
      mu = gene_rate * exp(cell_offsets)
      extra_signal = rnbinom(length(cells_toAdd_signal), mu = mu, size = 1 / gene_dispersion)
      
      # write the final, inflated expression into the sampled cells
      counts_inflated[gene_sample, cells_toAdd_signal] = 1 + extra_signal
    }
  }
  
  return(list(counts_inflated = counts_inflated, cell_info = cell_info))
}

# ============================================================================
# Given spatial coordinates, pairs each Sender ("ligand") spot with its
# nearest available Receiver ("receptor") spot, subject to a minimum distance
# threshold, then re-filters those pairs using their real pairwise distance.
#
# Arguments:
#   spatial_coords              - dataframe with 2 columns of coordinates, rownames = cell names
#   ligand_spots                - candidate Sender cell/spot names
#   receptor_spots               - candidate Receiver cell/spot names
#   distance_threshold_dataset   - initial threshold: only neighbors above this distance are considered
#   distance_post_filtering_lower/upper - after computing the real neighbor-neighbor
#                                          distance, cells outside this [lower, upper]
#                                          range are dropped
#
# Returns a list with:
#   metadata   - cell metadata with a Celltype column set to Sender / Receiver / Other
#   neighbors  - dataframe mapping each Sender cell to its paired Receiver cell
#                (the number of Sender and Receiver cells is always equal)
#   distance_Sender_toClosest_Receiver - average distance from each Sender to its Receiver
# ============================================================================
find_neighboring_spots = function(spatial_coords, ligand_spots, receptor_spots, distance_threshold_dataset,
                                  distance_post_filtering_lower, distance_post_filtering_upper)
{
  # ------------------------------------------------------------------
  # Step 1: pairwise distance matrix restricted to (Sender columns x Receiver rows)
  # ------------------------------------------------------------------
  df = st_as_sf(spatial_coords, coords = 1:2)
  dm = st_distance(df)
  
  colnames(dm) = rownames(spatial_coords)
  rownames(dm) = rownames(spatial_coords)
  
  # ligand_spots and receptor_spots never overlap, so this cleanly restricts
  # the matrix to Sender columns and Receiver rows
  dm = dm[, ligand_spots] %>% as.data.frame
  dm = dm[receptor_spots, ] %>% as.data.frame
  
  # ------------------------------------------------------------------
  # Step 2: for each Sender, greedily assign the closest still-unclaimed
  # Receiver that is above distance_threshold_dataset away
  # ------------------------------------------------------------------
  neighbors           = vector(length = ncol(dm))
  neighbors_dist       = vector(length = ncol(dm))
  receiver_cell_names = rownames(dm)
  
  for (index in seq_len(ncol(dm)))
  {
    # 1. distances from this Sender ('index') to every candidate Receiver
    dist_vec = dm[, index]
    
    # 2. mask out anything closer than the threshold, so it can't be picked below
    dist_vec[dist_vec < distance_threshold_dataset] = Inf
    
    # 3. take the closest in-range Receiver (NA/Inf if none is in range)
    c    = receiver_cell_names[which.min(dist_vec)]
    dist = dist_vec[which.min(dist_vec)]
    
    # Make sure different Sender cells don't "see" the same Receiver cell:
    # if this Receiver is already taken, walk down the sorted distance list
    # until we find one that isn't. This is a no-op when index == 1.
    n = 1
    while (any(c %in% neighbors))
    {
      c    = receiver_cell_names[order(dist_vec)[n]]
      dist = dist_vec[order(dist_vec)[n]]
      
      n = n + 1
      if (n > 5000) { message("while loop in fnc find_neighboring_spots taking too long"); break }
    }
    
    neighbors[index]      = c
    neighbors_dist[index] = dist
  }
  
  neighbors = as.list(neighbors)
  names(neighbors)      = colnames(dm) 
  names(neighbors_dist) = colnames(dm)
  
  # ------------------------------------------------------------------
  # Step 3: the greedy pairing above doesn't reflect the *real* distance from
  # each Sender to its closest Receiver, so recompute the real distances and
  # filter pairs down to a manually chosen [lower, upper] range
  # ------------------------------------------------------------------
  coords_sender   = metadata[unlist(neighbors) %>% names, c("x", "y")]
  coords_receiver = metadata[unlist(neighbors), c("x", "y")]
  
  dist_matrix = proxy::dist(coords_sender, coords_receiver, method = "Euclidean") %>% as.matrix
  
  # for each Sender, find its closest Receiver according to the real distance
  avg_dist_k1 = apply(dist_matrix, 1, function(row) {
    list(distance = sort(row)[1], closest_Receiver = names(sort(row)[1]))
  })
  
  # keep only pairs whose real distance falls inside [lower, upper]
  in_range    = lapply(avg_dist_k1, function(x) x$distance > distance_post_filtering_lower & x$distance < distance_post_filtering_upper) %>% unlist
  avg_dist_k1 = avg_dist_k1[in_range]
  
  neighbors = lapply(avg_dist_k1, function(x) x$closest_Receiver) %>% as.data.frame
  
  # some Receivers end up duplicated (expected, since we recomputed the closest
  # neighbor above) - keep only the first occurrence of each
  keep        = which(!duplicated(neighbors %>% unlist))
  neighbors   = neighbors[keep]
  avg_dist_k1 = avg_dist_k1[keep]
  
  neighbors_dist = lapply(avg_dist_k1, function(x) x$distance) %>% unlist
  
  # ------------------------------------------------------------------
  # Step 4: write the final Sender / Receiver / Other labels back to metadata
  # ------------------------------------------------------------------
  metadata %<>% mutate(Celltype = ifelse(Cell_ID %in% names(avg_dist_k1), "Sender", "Other"))
  metadata$Celltype[which(metadata$Cell_ID %in% unlist(neighbors))] = "Receiver"
  
  return(list(metadata = metadata, neighbors = neighbors,
              distance_Sender_toClosest_Receiver = unlist(neighbors_dist)))
}


# ============================================================================
# Estimates per-gene mean and dispersion using edgeR, used as the basis for
# the negative-binomial expression sampling in semi_simulate().
# ============================================================================
estimate_params_edgeR = function(counts, metadata, mm)
{
  dge = DGEList(counts = counts, samples = metadata)
  
  # estimate per-gene dispersion
  dge = estimateDisp(dge, design = mm)
  dge = edgeR::calcNormFactors(dge)
  
  # estimate mean expression per celltype group
  offset = edgeR::getOffset(dge)
  logmeans = edgeR::mglmOneWay(dge$counts, offset = offset, design = mm,
                               dispersion = dge$tagwise.dispersion)
  
  means_perCT = exp(logmeans$coefficients) %>% as.data.frame()
  colnames(means_perCT) = colnames(mm) %>% gsub("Celltype", "", .)
  means_perCT$gene_names = rownames(means_perCT)
  
  return(list(dge = dge, means_perCT = means_perCT, offset = offset))
}



# ============================================================================
# Generates an AveLogCPM diagnostic plot (original vs. inflated counts) for an
# explicit list of PCE_Sender/PCE_Receiver combinations, to visually
# confirm how much signal the semi-simulation added. 
#
# Arguments:
#   counts                - original (pre-inflation) count matrix
#   master_lst            - master_lst_diagnosticPlots (per-parameter-
#                            combination Sender/Receiver genes + precomputed
#                            inflated_counts_aveLogCPM)
#   metadata               - a representative cellmetadata dataframe (any one
#                            parameter combination's is fine - Sender/Receiver
#                            base Celltype assignment doesn't vary across
#                            combinations for a given dataset/strategy)
#   PCE_Sender_targets      - which PCE_Sender values to diagnose (a small,
#                            deliberately curated list, not the full config
#                            grid - see the call site for how these are chosen)
#   PCE_Receiver_targets    - same, for PCE_Receiver
#   dataset, CT_toPlot     - as before
#
# One plot is produced per (PCE_Sender_targets x PCE_Receiver_targets)
# combination
# ============================================================================
compute_diagnostic_plots = function(counts, master_lst, metadata, PCE_Sender_targets, PCE_Receiver_targets, dataset, CT_toPlot)
{
  plot_avelogcpm_fixed_PCE = list()
  
  # ---------------------------------------------------------------------
  # 1. Baseline: original AveLogCPM for the target cell types
  # ---------------------------------------------------------------------
  target_cells = metadata %>% filter(Celltype %in% CT_toPlot) %>% pull(Cell_ID)
  original_counts_subset = counts[, target_cells]
  original_aveLogCPM = aveLogCPM(original_counts_subset)
  
  # ---------------------------------------------------------------------
  # 2. One plot per (PCE_Sender, PCE_Receiver) combination
  # ---------------------------------------------------------------------
  PCE_grid = expand.grid(PCE_Sender = PCE_Sender_targets, PCE_Receiver = PCE_Receiver_targets)
  
  for (row in 1:nrow(PCE_grid))
  {
    PCE_Sender_target   = PCE_grid$PCE_Sender[row]
    PCE_Receiver_target = PCE_grid$PCE_Receiver[row]
    
    # master_lst is keyed by names like "PCE_Sender_X_PCE_Receiver_Y_indexLR_Z" -
    # take any one matching indexLR entry for this (PCE_Sender, PCE_Receiver) pair
    matches = grep(paste0("^PCE_Sender_", PCE_Sender_target, "_PCE_Receiver_", PCE_Receiver_target, "_indexLR_"), names(master_lst))
    if (length(matches) == 0) next
    param_key = names(master_lst)[matches[1]]
    
    # genes to highlight as L/R (also covers multi-subunit cases like R1_R2)
    LR_genes_color = c(master_lst[[param_key]]$Sender, master_lst[[param_key]]$Receiver) %>% str_split(., "_") %>% unlist
    
    # AveLogCPM of the inflated counts (precomputed and stored to save memory,
    # since we don't want to carry the entire inflated count matrix around)
    inflated_aveLogCPM = master_lst[[param_key]]$inflated_counts_aveLogCPM
    
    df_plot = data.frame(original_counts = original_aveLogCPM,
                         avelogcpm = inflated_aveLogCPM,
                         is_LR = names(inflated_aveLogCPM) %in% LR_genes_color)
    
    plot_avelogcpm_fixed_PCE[[param_key]] =
      ggplot(df_plot, aes(x = original_counts, y = avelogcpm, color = is_LR)) +
      geom_point(size = 0.5) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
      ggtitle(dataset, subtitle = paste0("PCE_S: ", PCE_Sender_target, " | PCE_R: ", PCE_Receiver_target)) +
      xlab("aveLogCPM (Original)") +
      ylab("aveLogCPM (Inflated)") +
      theme_light() +
      theme(
        plot.title  = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12)
      )
  }
  
  return(list(avelogcpm = plot_avelogcpm_fixed_PCE))
}

# ============================================================================
# Converts a JSON-derived nested list (as read by jsonlite::read_json) back
# into dataframes/vectors, handling a few fields specially.
# ============================================================================
convert_json_to_df = function(json_lst)
{
  out = lapply(names(json_lst), function(name)
  {
    lst = json_lst[name]
    
    if (name == "metadata") {
      # rebuild the metadata dataframe, with Cell_ID as rownames and x/y as numeric
      x = lst[[1]] %>% do.call(rbind.data.frame, .)
      rownames(x) = x$Cell_ID
      x = x %>% dplyr::mutate(across(c(x, y), as.numeric))
      
    } else if (name == "neighbor_cells") {
      x = lst$neighbor_cells %>% as.list %>% do.call(cbind.data.frame, .)
      
    } else if (name == "average_percentageCells_expressingLR") {
      x = lst %>% unlist
      
    } else if (name == "average_distance_SenderReceiver") {
      x = lst %>% unlist
    }
    
    return(x)
  }) %>% setNames(., names(json_lst))
  
  return(out)
}