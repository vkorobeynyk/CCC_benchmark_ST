# Semi-simulation framework
# one can semi-simulate several combinations of CT-CT pairs
semi_simulate = function(counts, simulated_interactions_lst ,genemetadata,  metadata , combination_CT, FC, n_neighbors, df_neighbors)
{
  n_neighbors_lst = list()
  counts_inflated = counts
  cell_info = list()
  means_perCT = genemetadata$mean
  
  for(comb_CT in combination_CT)
  {
    tmp_var1 = str_split(comb_CT,"_")[[1]]
    CT1 = tmp_var1[1]
    CT2 = tmp_var1[2]
    
    L_sample = simulated_interactions_lst[[comb_CT]] %>% str_split("_") %>% lapply(.,"[[",1) %>% as.character %>% setdiff("subunit") # remove subunit string from the L and R vectors
    R_sample = simulated_interactions_lst[[comb_CT]] %>% str_split("_") %>% lapply(.,"[[",2) %>% as.character %>% setdiff("subunit")
    
    # iterate over cell type
    for(tmp_CT in c("CTsender","CTreceiver"))  
    {
      if (tmp_CT == "CTsender" ) { # If celltype is sender -> add signal to all cells
        genes_to_sample = L_sample
        CT = CT1
        cells_toAdd_signal = colnames(counts)[which(metadata$Celltype == CT)]
      } else if (tmp_CT == "CTreceiver") {
        genes_to_sample = R_sample
        CT = CT2
        cells_toAdd_signal = colnames(counts)[which(metadata$Celltype == CT)] # CT2 receiver cells that we will add signal to were already precomputed with find_neighboring_spots, so we add signal to all cells
        }
      
      
      #save the inflated cells for estimating mean
      cell_info[[comb_CT]][["cells_signalAdded"]][[CT]] = cells_toAdd_signal
      
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
  } 
  
  return(list(counts_inflated = counts_inflated , cell_info = cell_info))
}

# Based on spatial coordinates data, this function selects neighboring cells 
find_neighboring_spots = function(metadata , spatial_coords, n_neighbors, ligand_spots, receptor_spots, remove_spots)
{
  # n_neighbors -> amount of neighboring spots/cells
  # ligand_spots -> spots/cells around which to select neighbors 
  # spatial_coords -> spatial coordinates. dataframe with 1 and 2 column being the coordinates
  # metadata -> dataframe with optional metadata
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
compute_diagnostic_plots = function(counts , master_lst, FC_param, n_neighbors_param, dataset, metadata , CT_toPlot)
{
  plot_avelogcpm_fixed_n_neighbors = list()
  
  # selects cells belonging to the celltype indicated by CT_toPlot
  filtered_raw_metadata =  filter(metadata, Celltype %in% CT_toPlot) 
  filtered_original_counts = counts[,filtered_raw_metadata$Cell_ID]
  filtered_original_counts_aveLogCPM = aveLogCPM(filtered_original_counts)
  
  ###############################################
  ### Plot AveLogCPM having n_neighbors fixed ###
  ###############################################
  for(n_neighbors in n_neighbors_param)
  {
    n = grep(paste0("^",n_neighbors,"$"), names(master_lst) %>% str_split("_") %>% lapply("[[", 5))
    # only select those files with correct FC
    n = intersect(n, which((names(master_lst) %>% str_split("_") %>% lapply("[[", 2)) %in% FC_param)) 
    # Get index of file with lowest FC
    min_FC_index = n[str_split( names(master_lst[n]),"_") %>% lapply(., "[[", 2) %>% which.min]
    max_FC_index = n[str_split( names(master_lst)[n],"_") %>% lapply(., "[[", 2) %>% which.max]
    
    ######################## Plot change in expression magnitude
    for(x in names(master_lst[n]))
    {
      # Select what genes to plot
      genes_to_plot = vector()
      for(i in CT_toPlot) {genes_to_plot = append(genes_to_plot , master_lst[[x]][[i]])}
      genes_to_plot = genes_to_plot %>% unlist %>% unique
      
      inflated_counts_aveLogCPM = master_lst[[x]][["inflated_counts_aveLogCPM"]]
      
      current_FC = str_split( x,"_") %>% lapply(., "[[", 2) %>% unlist
      df = data.frame(x = filtered_original_counts_aveLogCPM , y = inflated_counts_aveLogCPM, is_LR =  names(inflated_counts_aveLogCPM) %in% genes_to_plot)
      plot = ggplot(df,aes(x = x , y = y , color = is_LR)) + 
        geom_point(size = 0.5) + 
        ggtitle(paste0("n_neighbors = " ,n_neighbors , " dataset = ",dataset, " CT = ",paste(CT_toPlot, collapse = " "))) +
        xlab("aveLogCPM original counts") +
        ylab(paste("aveLogCPM FC=",current_FC)) +
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
plot_FCafter_semisimulation = function(vec, theoreticalFC, n_neighbors)
{
  theoreticalFC = theoreticalFC %>% unname()
  df = data.frame(gene = names(vec), value = vec)
  median = median(df$value) %>% round(.,2)
  plot = ggplot(df, aes(x = gene , y = value) ) + 
    geom_point() +
    geom_hline(yintercept=theoreticalFC, linetype="dashed", color = "red", linewidth = 1)  + 
    geom_hline(yintercept=median, linetype="dashed", color = "blue", linewidth = 1)  + 
    ggtitle(paste0("n_neighbors=",n_neighbors , " | theoretical FC=",theoreticalFC , " | real FC median=",median)) +
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
