# Semi-simulation framework
# one can semi-simulate several combinations of CT-CT pairs
semi_simulate = function(counts, simulated_interactions_lst ,fraction_cells_expressingR,fraction_cells_expressingL , genemetadata,  metadata , FC, df_neighbors)
{
  counts_inflated = counts
  cell_info = list()
  means_perCT = genemetadata$mean
  
  tmp_var1 = str_split("CT1_CT2","_")[[1]]
  CT1 = tmp_var1[1]
  CT2 = tmp_var1[2]
  
  L_sample = simulated_interactions_lst$CT1_CT2$ligand %>% str_split("_") %>% unlist
  R_sample = simulated_interactions_lst$CT1_CT2$receptor %>% str_split("_") %>% unlist
  
  # careful with the order of CTsender and CTreceiver
  # as I am using same variable cells_toAdd_signal, receiver has to come after sender
  for(tmp_CT in c("CTsender","CTreceiver"))  
  {
    # We add signal to either sender cells or receiver cells according to N_cells_expressing argument
    
    if (tmp_CT == "CTsender" ) { 
      genes_to_sample = L_sample
      CT = CT1
      N_cells_fromFraction = (length(colnames(df_neighbors)) * fraction_cells_expressingL) %>% floor
      cells_toAdd_signal = sample(colnames(df_neighbors), N_cells_fromFraction)
      
    } else if (tmp_CT == "CTreceiver") {
      genes_to_sample = R_sample
      CT = CT2
      
      # sample receiver cells according to ratio of fraction_cells_expressingL / fraction_cells_expressingR
      if(fraction_cells_expressingR/fraction_cells_expressingL == 1)
      {
        cells_toAdd_signal = df_neighbors[,cells_toAdd_signal] %>% unlist %>% unname
      } else if(fraction_cells_expressingR/fraction_cells_expressingL > 1) # in this case we have more Receiver cells than senders
      {
        n = (N_cells_fromFraction * fraction_cells_expressingR/fraction_cells_expressingL) %>% floor
        
        # select the real neighbors
        cells_toAdd_signal_realNeighbors = df_neighbors[,cells_toAdd_signal] %>% 
          na.omit %>%
          t %>% 
          as.data.frame %>% 
          sample_n(N_cells_fromFraction)
        # select cells which are not neighbors
        cells_toAdd_signal_nonrealNeighbors = df_neighbors[,which(!colnames(df_neighbors) %in% cells_toAdd_signal)]  %>% 
          na.omit %>%
          t %>% 
          as.data.frame %>% 
          sample_n(n-N_cells_fromFraction)
        
        cells_toAdd_signal = rbind(cells_toAdd_signal_realNeighbors, cells_toAdd_signal_nonrealNeighbors) %>% unlist %>% unname
      } else if(fraction_cells_expressingR/fraction_cells_expressingL < 1) # Sender cells more than receiver
      {
        cells_toAdd_signal = df_neighbors[,cells_toAdd_signal] %>% 
          na.omit %>%
          t %>% 
          as.data.frame %>% 
          sample_n((N_cells_fromFraction * fraction_cells_expressingR/fraction_cells_expressingL) %>% floor) %>% 
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
find_neighboring_spots = function(spatial_coords, max_N_neighbors, ligand_spots, receptor_spots, remove_spots)
{
  # max_N_neighbors -> amount of neighboring spots/cells
  # ligand_spots -> spots/cells around which to select neighbors 
  # spatial_coords -> spatial coordinates. dataframe with 1 and 2 column being the coordinates, rownames must be cellnames
  # remove_spots -> logical if to remove spots based on mean distance
  # CT1 -> name of celltype sender
  # CT2 -> name of celltype receiver
  
  # This function returns a metadata dataframe with a Celltype colum specifying Ligand, Receptor or Other spots
  
  # compute distances between all possible cells/spots
  df = st_as_sf(spatial_coords, coords=1:2)
  dm = st_distance(df)
  
  colnames(dm) = rownames(spatial_coords)
  rownames(dm) = rownames(spatial_coords)
  
  # select columns belonging to only ligand spots
  # ligand_spots and receptor_spots have no common entries
  dm = dm[,ligand_spots] %>% as.data.frame
  dm = dm[receptor_spots,] %>% as.data.frame
  
  # Select receptors that are closest to a ligand
  # get distance between ligand spot and neighboring receptor spot
  neighbors = vector(length = ncol(dm))
  neighbors_dist = vector(length = ncol(dm))
  receiver_cell_names = rownames(dm)
  for(index in seq_along(1:ncol(dm)))
  {
    c = receiver_cell_names[dm[,index] %>% order %>% extract(1:max_N_neighbors)]
    dist = dm[,index] %>% sort %>% extract(1:max_N_neighbors)
    
    # Make sure that different ligands dont "see" the same receptor
    n = 1
    while(c %in% neighbors)
    {
      c = receiver_cell_names[dm[,index] %>% order %>% extract(n)]
      dist = dm[,index] %>% sort %>% extract(n)
      n = n+1
      if(n > 5000) {message("while loop in fnc find_neighboring_spots taking too long") ; break}
    }
    neighbors[index] = c
    neighbors_dist[index] = dist
  }
  neighbors = as.list(neighbors)
  names(neighbors) = colnames(dm) # set names to ligand
  names(neighbors_dist) = colnames(dm) # set names to ligand
  
  ####################################################################
  ### Remove Receiver cells that are "far" away from Sender cells  ###
  if(remove_spots)
  {
    min_dist = neighbors_dist %>% unlist %>% mean()
    neighbors_dist = neighbors_dist[neighbors_dist < min_dist]
    neighbors = neighbors[names(neighbors) %in% names(neighbors_dist)]
    
    ligand_spots = names(neighbors)
    receiver_spots = unlist(neighbors)
  }
  
  # Add Celltype information to metadata file
  metadata %<>% mutate(Celltype = ifelse(Cell_ID %in% ligand_spots , "CT1", "Other"))
  metadata$Celltype[which(metadata$Cell_ID %in% receiver_spots)] = "CT2"
  
  return(list(metadata = metadata, neighbors = neighbors %>% as.data.frame))
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
    
    # as max_N_neighbors changes the amount of cells that I inflated counts into, I have to comptue avelogCPM for every different n_neighbor param
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

# As theoretical FC that we apply in the semi-simulation actually doesnt represent the practical FC that the data will be transformed with, generate a
# plot with real FC after semi-simulation
plot_FCafter_semisimulation = function(vec, indexLR,theoreticalFC, max_N_neighbors)
{
  indexLR = indexLR %>% unname()
  df = data.frame(gene = names(vec), value = vec)
  median = median(df$value) %>% round(.,2)
  plot = ggplot(df, aes(x = gene , y = value) ) + 
    geom_point() +
    geom_hline(yintercept=theoreticalFC, linetype="dashed", color = "red", linewidth = 1)  + 
    geom_hline(yintercept=median, linetype="dashed", color = "blue", linewidth = 1)  + 
    ggtitle(paste0("max_N_neighbors=",max_N_neighbors , " | indexLR=",indexLR , " | theoreticalFC=", theoreticalFC, " | real FC median=",median)) +
    xlab("LR index") +
    ylab("FC after simulation") +
    theme(axis.text.x=element_blank(), #remove x axis labels
          plot.title = element_text(size=8)  , 
          axis.text.y = element_text(size = 8)
    )
  return(list(plot = plot,theoreticalFC = theoreticalFC,  max_N_neighbors = max_N_neighbors,FC_real_median = median))
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
    } else if(name == "average_percentageCells_expressingLR") {x = lst %>% unlist}
    return(x)
  }) %>% setNames(., names(json_lst))
  return(out)
}

cumulative_expression = function(gene, seed,counts, cellmetadata) 
{
  set.seed(seed)
  
  # select randomly same number of cells for every celltype
  metadata = cellmetadata$metadata %>%
    group_by(Celltype) %>%
    sample_n(size = 50) %>%
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
  df_long$expr_sum =  df_long$expr_to  %>% as.numeric# for RECEPTORS
  #df_long$expr_sum =  df_long$expr_from # for LIGANDS
  #df_long$expr_sum =  df_long$expr_from + df_long$expr_to # for RECEPTORS + LIGANDS
  
  celltypes = unique(metadata$Celltype)
  
  # you can change bin_size if needed
  bin_size = 50
  max_dist = 500 # max(df_long$distance)
  
  dist_bins = seq(0, max_dist, by = bin_size)
  
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
