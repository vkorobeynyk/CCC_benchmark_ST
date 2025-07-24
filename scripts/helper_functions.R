# Semi-simulation framework
# one can semi-simulate several combinations of CT-CT pairs
semi_simulate = function(counts, simulated_interactions_lst ,N_cells_expressingR,N_cells_expressingL , genemetadata,  metadata , combination_CT, FC, df_neighbors)
{
  n_neighbors_lst = list()
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
      cells_toAdd_signal = sample(colnames(df_neighbors), N_cells_expressingL)
      
    } else if (tmp_CT == "CTreceiver") {
      genes_to_sample = R_sample
      CT = CT2
      # select receiver cells rowwise so that every sender cells is surrounded by same (+-1) receiver cell
      cells_toAdd_signal = df_neighbors[,cells_toAdd_signal] %>% 
        t %>% 
        as.vector %>% 
        extract(1:N_cells_expressingR)
      }
    
    
    #save the inflated cells for estimating mean
    cell_info[["CT1_CT2"]][["cells_signalAdded"]][[CT]] = cells_toAdd_signal
    
    # Iterate over every gene (L/R) depending on the CT and inflate expression
    for(gene_sample in genes_to_sample)
    {
      # set all the expression for this celltype to 0
      counts_inflated[gene_sample ,cells_toAdd_signal] = 0
      
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
find_neighboring_spots = function(spatial_coords, n_neighbors, ligand_spots, receptor_spots, remove_spots)
{
  # n_neighbors -> amount of neighboring spots/cells
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
  dm = dm[,ligand_spots] %>% as.data.frame
  dm = dm[receptor_spots,] %>% as.data.frame
  
  # select neighbor spots
  neighbors = apply(dm, 2, function(x, cell_names = rownames(dm)) {
    n = x %>% order %>% extract(2:(n_neighbors+1))
    n = cell_names[n]
  }) %>% as.data.frame %>% as.list
  
  # get distance between ligand spot and all neighboring receptor spots
  neighbors_dist = apply(dm, 2, function(x, cell_names = rownames(dm)) {
    n = x %>% sort %>% extract(2:(n_neighbors+1))
  }) %>% as.data.frame %>% as.list
  
  ####################################################################
  ### Remove Receiver cells that are "far" away from Sender cells  ###
  if(remove_spots)
  {
    min_dist = neighbors_dist %>% unlist %>% mean()
    for(ligand in ligand_spots)
    {
      # remove receiver cells that are above minimum distance
      x = neighbors_dist[[ligand]] < min_dist
      neighbors_dist[[ligand]] = neighbors_dist[[ligand]][x]
      neighbors[[ligand]] = neighbors[[ligand]][x]
      
      # as we want to have n_neighbors receiver cells around sender, we need to make sure that our receiver cells are not senders.
      # here remove those receiver cells that are actually also senders
      x = !neighbors[[ligand]] %in% ligand_spots
      neighbors[[ligand]] = neighbors[[ligand]][x]
      neighbors_dist[[ligand]] = neighbors_dist[[ligand]][x]
      
      # if sender has less neighbors than n_neighbors, remove
      if(length(neighbors_dist[[ligand]]) < n_neighbors)
      {
        neighbors_dist[[ligand]] = NULL
        neighbors[[ligand]] = NULL
        ligand_spots = ligand_spots %>% setdiff(ligand)
      }
    }
  }
  
  # Add Celltype information to metadata file
  metadata %<>% mutate(Celltype = ifelse(Cell_ID %in% ligand_spots , "CT1", "Other"))
  metadata$Celltype[metadata$Cell_ID %in% unlist(unname(neighbors)) & metadata$Celltype != "CT1"] = "CT2"
  
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
compute_diagnostic_plots = function(counts , master_lst, n_neighbors_param, dataset , CT_toPlot)
{
  plot_avelogcpm_fixed_n_neighbors = list()
  
  ###############################################
  ### Plot AveLogCPM across n_neighbors param ###
  ###############################################
  for(n_neighbors in n_neighbors_param)
  {
    n = grep(paste0("^",n_neighbors,"$"), names(master_lst) %>% str_split("_") %>% lapply("[[", 5))
    # select the first. Here we are just selecting indexLR which doesnt matter for the plors
    n = n[1]
    
    ######################## Plot change in expression magnitude
    for(x in names(master_lst[n]))
    {
      # Select what genes to label as L_R
      LR_genes_color = c(master_lst[[x]]$CT1,master_lst[[x]]$CT2)
      
      # as n_neighbors changes the amount of cells that I inflated counts into, I have to comptue avelogCPM for every different n_neighbor param
      metadata = master_lst[[x]]$simulated_cellmetadata$metadata
      
      # selects cells belonging to the celltype indicated by CT_toPlot
      filtered_raw_metadata =  filter(metadata, Celltype %in% CT_toPlot) 
      filtered_original_counts = counts[,filtered_raw_metadata$Cell_ID]
      filtered_original_counts_aveLogCPM = aveLogCPM(filtered_original_counts)
      
      inflated_counts_aveLogCPM = master_lst[[x]]$inflated_counts_aveLogCPM
      
      df = data.frame(x = filtered_original_counts_aveLogCPM , y = inflated_counts_aveLogCPM, is_LR =  names(inflated_counts_aveLogCPM) %in% LR_genes_color)
      plot = ggplot(df,aes(x = x , y = y , color = is_LR)) + 
        geom_point(size = 0.5) + 
        ggtitle(paste0("n_neighbors = " ,n_neighbors , " dataset = ",dataset, " CT = ",paste(CT_toPlot, collapse = " "))) +
        xlab("aveLogCPM original counts") +
        ylab(paste("aveLogCPM")) +
        theme(plot.title = element_text(size=8) , 
              axis.text.x = element_text(size = 8) , 
              axis.text.y = element_text(size = 8))
      plot_avelogcpm_fixed_n_neighbors[[x]] = plot
    }
  }
  return(list(avelogcpm_fixed_n_neighbors = plot_avelogcpm_fixed_n_neighbors))
}

# As theoretical FC that we apply in the semi-simulation actually doesnt represent the practical FC that the data will be transformed with, generate a
# plot with real FC after semi-simulation
plot_FCafter_semisimulation = function(vec, indexLR,theoreticalFC, n_neighbors)
{
  indexLR = indexLR %>% unname()
  df = data.frame(gene = names(vec), value = vec)
  median = median(df$value) %>% round(.,2)
  plot = ggplot(df, aes(x = gene , y = value) ) + 
    geom_point() +
    geom_hline(yintercept=theoreticalFC, linetype="dashed", color = "red", linewidth = 1)  + 
    geom_hline(yintercept=median, linetype="dashed", color = "blue", linewidth = 1)  + 
    ggtitle(paste0("n_neighbors=",n_neighbors , " | indexLR=",indexLR , " | theoreticalFC=", theoreticalFC, " | real FC median=",median)) +
    xlab("LR index") +
    ylab("FC after simulation") +
    theme(axis.text.x=element_blank(), #remove x axis labels
          plot.title = element_text(size=8)  , 
          axis.text.y = element_text(size = 8)
    )
  return(list(plot = plot,theoreticalFC = theoreticalFC,  n_neighbors = n_neighbors,FC_real_median = median))
}

# converts json list to dataframes
# iterates over every index in list and converts it to df
# when the list only contains 1 integer, converts to vector
# 
convert_json_to_df = function(json_lst) 
{
  x = lapply(json_lst, function(x) 
  {
    x1 = lapply(x, unlist) %>% do.call(rbind.data.frame, .)
    if(length(x1) == 1) {x1 = x1 %>% unlist %>% unname}
    colnames(x1) = x[[1]] %>% names
    
    if(any(colnames(x1) %in% c("x","y"))) 
    {
      x1 = x1 %>% mutate(across(c(x,y), as.double))
    }
    return(x1)
  })
  return(x)
}
