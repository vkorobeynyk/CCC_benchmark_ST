general_ST_QCpipeline_R = function(counts, to_filter=FALSE)
{
  ##########################
  #### load data as sce ####
  ##########################
  
  sce = SingleCellExperiment(list(counts = counts))
  
  #########################
  #### QC b/ filtering ####
  #########################
  message("############# QC b/ filtering #############")
  #### Per cell and per feature QC
  
  # calculate MT genes again as some can be filtered in the previous step
  mito.genes <- grep(pattern = "^MT-", x = rownames(counts(sce)), value = TRUE)
  
  sce <- addPerCellQC(sce, 
                      subsets = list(mito = mito.genes))
  
  # plot counts vs features detected
  plt_sum_detected = plotColData(sce, x = "sum", y = "detected")
  print(plt_sum_detected)
  
  # Plot features with highest average expression across all cells, along with their expression in each individual cell
  #rowData(sce)$genes = rownames(sce)
  #print(plotHighestExprs(sce, exprs_values = "counts", n = 50, feature_names_to_plot = "genes"))
  
  #### Plot explanatory Variables
  
  sce <- logNormCounts(sce)
  
  #vars <- getVarianceExplained(sce)
  
  #print(plotExplanatoryVariables(vars, ylim = c(0, 200)))
  #### Filtering parameters
  
  libsize.drop_lower <- vector()
  libsize.drop_upper <- vector()
  feature.drop_lower <- vector()
  feature.drop_upper <- vector()
  mito.drop <- vector()
  
  # Detect Outliers
  libsize.drop_lower <- isOutlier(sce$total, nmads = 3, type = "lower", log = TRUE)
  libsize.drop_upper <- isOutlier(sce$total, nmads = 2, type = "higher", log = TRUE)
  
  feature.drop_lower <- isOutlier(sce$detected, nmads = 3, type = "lower", log = TRUE)
  feature.drop_upper <- isOutlier(sce$detected, nmads = 2, type = "higher", log = TRUE)
  
  mito.drop <- isOutlier(sce$subsets_mito_percent, nmads = 2.5, type = "higher")
  
  ### Quality metrics and outlier detection
  pass <- !(libsize.drop_lower | libsize.drop_upper | feature.drop_lower | feature.drop_upper | mito.drop)
  sce$qc_pass <- factor(ifelse(test = pass,
                               yes = 'QC_pass',
                               no = 'QC_fail'))
  
  #### Overall quality metrics before QC
  
  dt = colData(sce) %>% as.data.frame
  plt_hist_libsize_before_filtering = ggplot(dt, aes(x = total)) + 
    geom_histogram(bins = 50,fill="red",alpha=0.5) + 
    xlab("Library size") +
    ylab("Number of cells") +
    geom_vline(xintercept =  attr(libsize.drop_lower, 'thresholds')['lower'], col="blue")+
    geom_vline(xintercept =  attr(libsize.drop_upper, 'thresholds')['higher'], col="blue") +
    scale_x_log10()
  print(plt_hist_libsize_before_filtering)
  
  plt_hist_Nexpressedgenes_before_filtering = ggplot(dt, aes(x = detected)) + 
    geom_histogram(bins = 50,fill="red",alpha=0.5) + 
    xlab("Number of expressed genes") +
    ylab("Number of cells") +
    geom_vline(xintercept =  attr(feature.drop_lower, 'thresholds')['lower'], col="blue")+
    geom_vline(xintercept =  attr(feature.drop_upper, 'thresholds')['higher'], col="blue")+
    scale_x_log10()
  print(plt_hist_Nexpressedgenes_before_filtering)
  
  plt_hist_mitproportion_before_filtering = ggplot(dt, aes(x = subsets_mito_percent)) + 
    geom_histogram(bins = 50,fill="red",alpha=0.5) + 
    xlab("Mitochondrial proportion (%)") +
    ylab("Number of cells") +
    geom_vline(xintercept =  attr(mito.drop, 'thresholds')['higher'], col="blue")+
    scale_x_log10()
  print(plt_hist_mitproportion_before_filtering)
  
  #### QC-pass and QC-fail
  
  print(plotColData(sce, x = "sum", y = "detected", colour_by = 'qc_pass') +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position="right", aspect.ratio = 1))
  
  #### QC-pass only
  
  print(plotColData(sce[,colData(sce)$qc_pass == "QC_pass"],
                    x = "sum", y = "detected", colour_by = 'qc_pass') +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position="right", aspect.ratio = 1))
  
  #### QC-fail only
  
  print(plotColData(sce[,colData(sce)$qc_pass == "QC_fail"],
                    x = "sum", y = "detected",colour_by = 'qc_pass') +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position="right", aspect.ratio = 1))
  
  Sys.sleep(2) # added here because the previous plot appears after the QC a/filtering messaage
  #########################
  #### QC a/ filtering ####
  #########################
  message("############# QC a/ filtering #############")
  
  if(to_filter) {sce_filtered = sce[,pass]} else {sce_filtered = sce}
  
  # Plot features with highest average expression across all cells, along with their expression in each individual cell
  #print(plotHighestExprs(sce_filtered, exprs_values = "counts", n = 50, feature_names_to_plot = "genes"))
  
  #### Overall quality metrics after QC
  
  dt = colData(sce_filtered) %>% as.data.frame
  plt_hist_libsize_after_filtering = ggplot(dt, aes(x = total)) + 
    geom_histogram(bins = 50,fill="red",alpha=0.5) + 
    xlab("Library size") +
    ylab("Number of cells")+
    scale_x_log10()
  print(plt_hist_libsize_after_filtering)
  
  plt_hist_Nexpressedgenes_after_filtering = ggplot(dt, aes(x = detected)) + 
    geom_histogram(bins = 50,fill="red",alpha=0.5) + 
    xlab("Number of expressed genes") +
    ylab("Number of cells")+
    scale_x_log10()
  print(plt_hist_Nexpressedgenes_after_filtering)
  
  plt_hist_mitproportion_after_filtering = ggplot(dt, aes(x = subsets_mito_percent)) + 
    geom_histogram(bins = 50,fill="red",alpha=0.5) + 
    xlab("Mitochondrial proportion (%)") +
    ylab("Number of cells")+
    scale_x_log10()
  print(plt_hist_mitproportion_after_filtering)
  
  Sys.sleep(2) # added here because the previous plot appears after the QC a/filtering messaage
  ###############################
  #### Additional QC metrics ####
  ###############################
  message("############# Additional QC metrics #############")
  
  plot(density(sce_filtered$total/1e3),
       xlab = "library size (thousands)",
       main = '')
  
  rug(sce_filtered$total/1e3)
  
  plot(y = sce_filtered$total/1e3,
       x = sce_filtered$subsets_mito_percent,
       pch = 20,
       ylab = 'library size (thousands)',
       xlab = 'mitochondrial proportion (%)',
       ### col = ac('black', 0.5)
  )
  
  plot(y = sce_filtered$total/1e3,
       x = sce_filtered$detected,
       pch = 20,
       ylab = 'library size (thousands)',
       xlab = 'number of genes',
       ### col = ac('black', 0.5)
  )
  return(list(sce_filtered = sce_filtered, 
              plt_sum_detected = plt_sum_detected,
              plt_hist_libsize_before_filtering = plt_hist_libsize_before_filtering,
              plt_hist_Nexpressedgenes_before_filtering = plt_hist_Nexpressedgenes_before_filtering,
              plt_hist_mitproportion_before_filtering = plt_hist_mitproportion_before_filtering,
              plt_hist_libsize_after_filtering = plt_hist_libsize_after_filtering,
              plt_hist_Nexpressedgenes_after_filtering = plt_hist_Nexpressedgenes_after_filtering,
              plt_hist_mitproportion_after_filtering = plt_hist_mitproportion_after_filtering))
}
