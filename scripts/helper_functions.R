# Semi-simulation framework
# one can semi-simulate several combinations of CT-CT pairs
semi_simulate = function(counts, simulated_interactions_lst ,fraction_cells_expressingR,fraction_cells_expressingL , genemetadata,  metadata , FC, df_neighbors)
{
  counts_inflated = counts
  cell_info = list()
  means_perCT = genemetadata$mean
  
  L_sample = simulated_interactions_lst$CT1_CT2$ligand %>% str_split("_") %>% unlist
  R_sample = simulated_interactions_lst$CT1_CT2$receptor %>% str_split("_") %>% unlist
  
  # careful with the order of CTsender and CTreceiver
  # as I am using same variable cells_toAdd_signal, receiver has to come after sender
  for(tmp_CT in c("CTsender","CTreceiver"))  
  {
    # We add signal to either sender cells or receiver cells according to N_cells_expressing argument
    
    if (tmp_CT == "CTsender" ) { 
      genes_to_sample = L_sample
      CT = "CT1"
      CT1_cells = (length(colnames(df_neighbors)) * fraction_cells_expressingL) %>% floor
      cells_toAdd_signal = sample(colnames(df_neighbors), CT1_cells)
      
    } else if (tmp_CT == "CTreceiver") {
      genes_to_sample = R_sample
      CT = "CT2"
      
      # sample receiver cells according to ratio of fraction_cells_expressingL / fraction_cells_expressingR
      if(fraction_cells_expressingR/fraction_cells_expressingL == 1)
      {
        cells_toAdd_signal = df_neighbors[,cells_toAdd_signal] %>% unlist %>% unname
      } else if(fraction_cells_expressingR/fraction_cells_expressingL > 1) # in this case we have more Receiver cells than senders
      {
        n = (length(colnames(df_neighbors)) * fraction_cells_expressingR) %>% floor
        
        # select the real neighbors
        cells_toAdd_signal_CT1_signal_Added_Neighbors = df_neighbors[,cells_toAdd_signal] %>% 
          na.omit %>%
          t %>% 
          as.data.frame %>% 
          sample_n(CT1_cells)
        # select CT2 cells
        cells_toAdd_signal_non_CT1_signal_Added_Neighbors = df_neighbors[,which(!colnames(df_neighbors) %in% cells_toAdd_signal)]  %>% 
          na.omit %>%
          t %>% 
          as.data.frame %>% 
          sample_n(n-CT1_cells)
        
        cells_toAdd_signal = rbind(cells_toAdd_signal_CT1_signal_Added_Neighbors, cells_toAdd_signal_non_CT1_signal_Added_Neighbors) %>% unlist %>% unname
      } else if(fraction_cells_expressingR/fraction_cells_expressingL < 1) # Sender cells more than receiver
      {
        cells_toAdd_signal = df_neighbors[,cells_toAdd_signal] %>% 
          na.omit %>%
          t %>% 
          as.data.frame %>% 
          sample_n((length(colnames(df_neighbors)) * fraction_cells_expressingR) %>% floor) %>% 
          unlist %>% unname
      }
    }
    
    #save the inflated cells for estimating mean
    cell_info[["CT1_CT2"]][["cells_signalAdded"]][[CT]] = cells_toAdd_signal
    
    # Iterate over every gene (L/R) depending on the CT and inflate expression
    for(gene_sample in genes_to_sample)
    {
      # set all the expression for this celltype (CT1 or CT2) to 0
      counts_inflated[gene_sample ,metadata$Celltype %>% grepl(CT,.)] = 0
      
      gene_mean = means_perCT[grep(paste("^",gene_sample,"$", sep=""),  means_perCT$gene_names),] %>% .[CT] %>% as.numeric()
      
      gene_dispersion = genemetadata$disp %>% subset(gene == gene_sample) %>% select(edgeR_dispersion) %>% as.numeric
      mu = gene_mean * FC
      x1 = rnbinom(1000, mu = mu, size = 1/gene_dispersion) # shape parameter of the gamma mixing distribution
      # replace the expression for the CT according to the sampled values
      # in case there are not enough sampled values > 0 then sample from the > 0 values with replacement
      if(length(x1[x1>0]) > length(cells_toAdd_signal)) {x1 = sample(x1[x1>0] , length(cells_toAdd_signal), replace = F)
      } else if(all(x1 == 0)) {x1 = sample(1 , length(cells_toAdd_signal), replace = T)
      } else {x1 = sample(x1[x1>0] , length(cells_toAdd_signal), replace = T)}
      
      
      # Add the final expression to sampled zero cells
      counts_inflated[gene_sample ,cells_toAdd_signal] = x1
    }
  }
  
  return(list(counts_inflated = counts_inflated , cell_info = cell_info))
}

# Based on spatial coordinates data, this function selects neighboring cells 
find_neighboring_spots = function(spatial_coords, ligand_spots, receptor_spots, remove_spots, distance_threshold_dataset,
                                  distance_post_filtering_lower, distance_post_filtering_upper)
{
  # ligand_spots -> cells/spots around which to select neighbors 
  # spatial_coords -> spatial coordinates. dataframe with 1 and 2 column being the coordinates, rownames must be cellnames
  # remove_spots -> logical, if to remove spots based on mean distance, Recommended = TRUE
  # CT1 -> name of celltype sender
  # CT2 -> name of celltype receiver
  # distance_threshold_dataset -> initial threshold to find neighbors that are above this threshold
  # distance_post_filtering_lower -> after calculate the real neighbor-neighbor distance, this is the lower threshold for filtering
  # distance_post_filtering_upper -> after calculate the real neighbor-neighbor distance, this is the upper threshold for filtering
  
  # This function returns a metadata dataframe with a Celltype column specifying sender, receiver or Other cells and 
  # another neighbor_spots which is a dataframe with sender (names) and corresponding receiver cell
  # The amount of Sender/Receiver cells/spots is the same
  # average_distance_CT1_CT2_during_sampling is the average distance between CT1 and CT2 during phase of selecting the neighbor (based on distance_threshold_dataset param)
  # average_distance_CT1_CT2_after_sampling is the average distance for every CT1 to CT2
  
  # compute distances between all possible cells/spots
  df = st_as_sf(spatial_coords, coords=1:2)
  dm = st_distance(df)
  
  colnames(dm) = rownames(spatial_coords)
  rownames(dm) = rownames(spatial_coords)
  
  # select columns belonging to only ligand spots
  # ligand_spots and receptor_spots have no common entries
  dm = dm[,ligand_spots] %>% as.data.frame
  dm = dm[receptor_spots,] %>% as.data.frame
  
  # --------------------------------------------------------------- #
  # Select receivers that are closest to a sender
  # if the receiver has already sender pair, find the 2nds closest receiver etc...
  # the receiver must be above distance_threshold_dataset distance 
  # --------------------------------------------------------------- #
  neighbors = vector(length = ncol(dm))
  neighbors_dist = vector(length = ncol(dm))
  receiver_cell_names = rownames(dm)
  for(index in seq_along(1:ncol(dm)))
  {
    # 1. Extract the column of distances for the current 'index'
    dist_vec = dm[, index]
    
    # 2. Filter for distances within the range
    # We set everything else to Inf so they aren't picked by order/which.min
    #dist_vec[dist_vec < 50 | dist_vec > 250] = Inf
    dist_vec[dist_vec < distance_threshold_dataset] = Inf
    
    # 3. Extract the name of the closest cell in that range
    # If no cell is in range, 'c' will be NA or point to an Inf value
    c = receiver_cell_names[which.min(dist_vec)]
    
    dist = dist_vec[which.min(dist_vec)]
    
    # Make sure that different ligands dont "see" the same receptor
    # This doesnt do anything when index == 1
    n = 1
    while(any(c %in% neighbors))
    {
      # order cells according to distance
      # iterate over the ordered list to find receiver that has no sender associated to it
      c = receiver_cell_names[order(dist_vec)[n]]
      dist = dist_vec[order(dist_vec)[n]]
      
      n = n+1
      if(n > 5000) {message("while loop in fnc find_neighboring_spots taking too long") ; break}
    }
    neighbors[index] = c
    neighbors_dist[index] = dist
  }
  neighbors = as.list(neighbors)
  names(neighbors) = colnames(dm) # set names to ligand
  names(neighbors_dist) = colnames(dm) # set names to ligand
  
  # --------------------------------------------------------------- #
  # Unfortunately the code above doesnt reflect the real distance from each sender to the closest receiver
  # Here I want to filter the neighbors found above based on the real distance 
  # Use upper and lower distance bound threshold manually selected
  # --------------------------------------------------------------- #
  ###### Compute average practical distance between all CT1 to all CT2 cells
  coords_ct1 = metadata[unlist(neighbors) %>% names, c("x", "y")]
  coords_ct2 = metadata[unlist(neighbors), c("x", "y")]
  
  dist_matrix = proxy::dist(coords_ct1, coords_ct2, method = "Euclidean")
  
  # Convert to a standard matrix if needed
  dist_matrix = as.matrix(dist_matrix)
  
  # Calculate the dist of the closest CT2 cells for each CT1 cell
  avg_dist_k1 = apply(dist_matrix, 1, function(row) {
    # Sort the distances and take the first
    return(list(distance = sort(row)[1],closest_CT2 = names(sort(row)[1])))
  })
  
  # filter based on distance
  avg_dist_k1 = avg_dist_k1[lapply(avg_dist_k1, function(x) {x$distance > distance_post_filtering_lower & x$distance < distance_post_filtering_upper}) %>% unlist]
  
  neighbors = lapply(avg_dist_k1, function(x) {x$closest_CT2}) %>% as.data.frame
  # some receiver are duplicated -> normal because we recomputed the closest neighbor, remove those
  n = which(!duplicated(neighbors %>% unlist))
  neighbors = neighbors[n]
  
  avg_dist_k1 = avg_dist_k1[n]
  neighbors_dist = lapply(avg_dist_k1, function(x) {x$distance}) %>% unlist
  
  # update metadata
  metadata %<>% mutate(Celltype = ifelse(Cell_ID %in% names(avg_dist_k1) , "CT1", "Other"))
  metadata$Celltype[which(metadata$Cell_ID %in% unlist(neighbors))] = "CT2"
  
  return(list(metadata = metadata, neighbors = neighbors, 
              distance_CT1_toClosest_CT2 = unlist(neighbors_dist)))
}

select_HighDensityRegion = function(metadata, neighbor_cells, dataset, sample_cells = 100)
{
  ##
  df = metadata
  
  df_binned = df %>%
    mutate(
      x_bin = cut(x, breaks = 4),
      y_bin = cut(y, breaks = 4)
    )
  
  bin_counts = df_binned %>%
    dplyr::count(x_bin, y_bin, Celltype)
  
  bin_wide = bin_counts %>%
    tidyr::pivot_wider(names_from = Celltype, values_from = n, values_fill = 0)
  
  bin_wide = bin_wide %>%
    mutate(co_score = CT1 * CT2)
  
  top_regions = bin_wide %>%
    dplyr::arrange(desc(co_score)) %>%
    mutate(cells = CT1 + CT2 + Other)
  
  #top_regions %>% arrange(desc(cells))
  
  if(dataset == "Visium_HD_HPC")
  {
    # from the high density region, select CT1 and CT2 cells according to sample_cells
    metadata$Celltype = "Other"
    ct1cells = metadata %>% filter(between(metadata$x, 3.5e+03, 4.87e+03) & between(metadata$y, 1.58e+04, 1.8e+04)) %>% sample_n(sample_cells) %>% pull(Cell_ID)
    metadata[ct1cells,"Celltype"] = "CT1"
    ct2cells = metadata %>% filter(between(metadata$x, 3.5e+03, 4.87e+03) & between(metadata$y, 1.58e+04, 1.8e+04) & Celltype != "CT1") %>% sample_n(sample_cells) %>% pull(Cell_ID)
    metadata[ct2cells,"Celltype"] = "CT2"
    
    # generate neighbor_cells dataframe without filtering for far away connetions
    neighbor_info = find_neighboring_spots(spatial_coords = metadata %>% select(c("x","y")), 
                                           ligand_spots = metadata %>% filter(Celltype == "CT1") %>% select(Cell_ID) %>% unlist %>% unname, 
                                           receptor_spots = metadata %>% filter(Celltype == "CT2") %>% select(Cell_ID) %>% unlist %>% unname,
                                           remove_spots = FALSE)
    
  } else {stop(paste0("You must specify coordinates for dataset ", dataset))}
  
  # make the amount of CT1 and CT2 cells equal 
  #x = table(metadata$Celltype)
  #if(x["CT1"] < x["CT2"])
  #{
  #  neighbor_cells = neighbor_cells[metadata$Cell_ID[which(metadata$Celltype == "CT1")]]
  #  metadata$Celltype = "Other"
  #  metadata[names(neighbor_cells), "Celltype"] = "CT1"
  #  metadata[neighbor_cells, "Celltype"] = "CT2"
  #} else {
  #  neighbor_cells = neighbor_cells[neighbor_cells %in% metadata$Cell_ID[which(metadata$Celltype == "CT2")]]
  #  metadata$Celltype = "Other"
  #  metadata[names(neighbor_cells), "Celltype"] = "CT1"
  #  metadata[unlist(neighbor_cells), "Celltype"] = "CT2"
  #}
  return(list(neighbor_cells = neighbor_info$neighbors,
              metadata = neighbor_info$metadata))
}


### Estimates mean and dispersion using edgeR
estimate_params_edgeR = function(counts , metadata, mm)
{
  dge = DGEList(counts = counts, samples = metadata)
  
  # Estimate disp
  dge = estimateDisp(dge , design = mm)
  dge = edgeR::calcNormFactors(dge)
  
  # estimating mu
  #centered.off = edgeR::getOffset(dge)  
  #centered.off = centered.off - mean(centered.off) 
  logmeans = edgeR::mglmOneWay(dge$counts, offset = 0, design = mm,
                               dispersion = dge$tagwise.dispersion) 
  
  means_perCT = exp(logmeans$coefficients) %>% as.data.frame()
  colnames(means_perCT) = colnames(mm)
  colnames(means_perCT) = gsub("Celltype","",colnames(means_perCT))
  means_perCT$gene_names = rownames(means_perCT)
  
  return(list(dge =  dge , means_perCT = means_perCT))
}

# Generate 2 plots:
# avelogcpm plot according to edgeR that shows how much signal we added to data
compute_diagnostic_plots = function(counts , master_lst, indexLR ,FC_nSenderCells , FC_nReceiverCells, dataset , CT_toPlot)
{
  plot_avelogcpm_fixed = list()
  
  # select some combinations for plotting
  params_grid = expand.grid(indexLR = indexLR, FC_nSenderCells = FC_nSenderCells, FC_nReceiverCells = FC_nReceiverCells) %>% sample_n(2)
  
  ########################################
  ### Plot AveLogCPM across grid_param ###
  ########################################
  for(i in 1:nrow(params_grid))
  {
    x = params_grid[i,] %>% as.numeric ;  names(x) = c("indexLR","FC_nSenderCells","FC_nReceiverCells")
    
    n = grep(paste0("FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"],"_indexLR_",x["indexLR"]), names(master_lst))
    file = names(master_lst)[n]
    # Select what genes to label as L_R
    LR_genes_color = c(master_lst[[file]]$CT1,master_lst[[file]]$CT2) %>% str_split(., "_") %>% unlist # for cases when we have R1_R2 subunits
    
    metadata = master_lst[[file]]$simulated_cellmetadata$metadata
    
    # selects cells belonging to the celltype indicated by CT_toPlot
    filtered_raw_metadata =  filter(metadata, Celltype %in% CT_toPlot) 
    filtered_original_counts = counts[,filtered_raw_metadata$Cell_ID]
    filtered_original_counts_aveLogCPM = aveLogCPM(filtered_original_counts)
    
    inflated_counts_aveLogCPM = master_lst[[file]]$inflated_counts_aveLogCPM
    
    df = data.frame(x = filtered_original_counts_aveLogCPM , y = inflated_counts_aveLogCPM, is_LR =  names(inflated_counts_aveLogCPM) %in% LR_genes_color)
    plot = ggplot(df,aes(x = x , y = y , color = is_LR)) + 
      geom_point(size = 2) + 
      ggtitle(paste0("indexLR_",x["indexLR"] , "_FC_nSenderCells_",x["FC_nSenderCells"], "_FC_nReceiverCells_", x["FC_nReceiverCells"] , " dataset = ",dataset, " CT = ",paste(CT_toPlot, collapse = " "))) +
      xlab("aveLogCPM original counts") +
      ylab(paste("aveLogCPM")) +
      theme(plot.title = element_text(size=8) , 
            axis.text.x = element_text(size = 8) , 
            axis.text.y = element_text(size = 8))
    plot_avelogcpm_fixed[[file]] = plot
  }
  return(list(avelogcpm = plot_avelogcpm_fixed))
}

# converts json list to dataframes
# iterates over every index in list and converts it to df
# when the list only contains 1 integer, converts to vector
convert_json_to_df = function(json_lst) 
{
  out = lapply(names(json_lst), function(name) 
  {
    lst = json_lst[name]
    # If second index of list is non existant, its a vector
    if(name == "metadata") { 
      x = lst[[1]] %>% do.call(rbind.data.frame, .)
      rownames(x) = x$Cell_ID
      x = x %>% dplyr::mutate(across(c(x,y), as.numeric))
    } else if(name == "neighbor_cells") {x = lst$neighbor_cells %>% as.list %>% do.call(cbind.data.frame,.)
    } else if(name == "average_percentageCells_expressingLR") {x = lst %>% unlist
    } else if(name == "N_Other_CT2_cells_expressing_ligand") {x = lst %>% unlist
    } else if(name == "N_Other_CT1_cells_expressing_receptor") {x = lst %>% unlist}
    return(x)
  }) %>% setNames(., names(json_lst))
  return(out)
}

cumulative_expression = function(gene, seed,counts, cellmetadata, test_receptors = TRUE ,N_cells = 50, distance_bin_size = 50, max_dist = 500) 
{
  set.seed(seed)
  
  # select randomly same number of cells for every celltype
  metadata = cellmetadata %>%
    group_by(Celltype) %>%
    sample_n(size = N_cells) %>%
    ungroup()
  
  counts_subset = counts[,metadata$Cell_ID]
  # compute distances between all possible cells/spots
  df = st_as_sf(metadata %>% select(x,y), coords=1:2)
  dm = st_distance(df)
  
  # named vector for direct lookup
  celltype_vec = setNames(metadata$Celltype, metadata$Cell_ID)
  
  # add celltype to cellnames for next step
  colnames(dm) = rownames(dm) = colnames(counts_subset) %>% stringr::str_c(., "-", celltype_vec)
  
  df_long = dm %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from") %>%
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "distance"
    ) %>% 
    dplyr::filter(distance != 0) # drop self-distance
  
  # add celltype column
  df_long$ctype_from = str_split(df_long$from, "-") %>% lapply(.,"[[",2) %>% unlist
  df_long$ctype_to   = str_split(df_long$to, "-") %>% lapply(.,"[[",2) %>% unlist
  df_long$from = str_split(df_long$from, "-") %>% lapply(.,"[[",1) %>% unlist # remove celltype info from cellnames
  df_long$to   = str_split(df_long$to, "-") %>% lapply(.,"[[",1) %>% unlist # remove celltype info from cellnames
  
  gene_expr = counts_subset[gene, ]
  
  df_long$expr_from = gene_expr[df_long$from] %>% as.numeric
  df_long$expr_to   = gene_expr[df_long$to] %>% as.numeric
  if(test_receptors) {df_long$expr_sum =  df_long$expr_to  %>% as.numeric # for RECEPTORS
  } else if(test_receptors == FALSE) {df_long$expr_sum =  df_long$expr_from %>% as.numeric  # for LIGANDS
  } else{message("specify to test either ligand or receptor gene"); break}
  
  
  
  #df_long$expr_sum =  df_long$expr_from # for LIGANDS
  #df_long$expr_sum =  df_long$expr_from + df_long$expr_to # for RECEPTORS + LIGANDS
  
  celltypes = unique(metadata$Celltype)
  
  dist_bins = seq(100, max_dist, by = distance_bin_size)
  
  curve_list = lapply(celltypes, function(sender) {
    
    ##############
    # FIX SENDER #
    ##############
    edges_ct = df_long %>% 
      filter(ctype_from == sender)
    
    ################
    # FIX RECEIVER #
    ################
    lapply(celltypes[celltypes != sender], function(receiver) {
      # Compute cumulative expression for each radius
      edges_ct2 = edges_ct %>% 
        filter(ctype_to == receiver)
      
      data.frame(
        sender_receiver = paste(sender, receiver),
        radius = dist_bins,
        total_expr = sapply(dist_bins, function(d) sum(edges_ct2$expr_sum[edges_ct2$distance <= d], na.rm = TRUE)
        )
      )
    })
  })
  
  distance_curve_ct = bind_rows(curve_list)
  
  p1 = ggplot(distance_curve_ct, aes(x = radius, y = total_expr, color = sender_receiver)) +
    geom_line(size = 1.2) +
    geom_point() +
    theme_classic(base_size = 14) +
    labs(
      x = "Distance",
      y = paste("Summed normalized expression of", gene),
      title = paste("Cumulative gene expression by radius:", gene)
    ) +
    scale_color_brewer(palette = "Dark2") 
  return(list(distance_curve_ct = distance_curve_ct, plot = p1))
}
