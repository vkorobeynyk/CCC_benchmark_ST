# ============================================================================
# general_ST_QCpipeline_R()
#
# Generic QC pipeline for spatial transcriptomics count data. Computes
# per-cell QC metrics (library size, number of detected genes, mitochondrial
# percentage), flags outlier cells using median-absolute-deviation (MAD)
# based thresholds, produces diagnostic QC plots before/after filtering, and
# optionally returns the filtered SingleCellExperiment.
#
# Arguments:
#   counts     - gene x cell count matrix
#   to_filter  - logical; if TRUE, cells failing QC are removed from the
#                returned sce_filtered object. If FALSE, all cells are kept
#                (QC metrics/plots are still computed either way).
#
# Returns a list with the filtered (or unfiltered) SingleCellExperiment and
# all diagnostic ggplot objects produced along the way.
# ============================================================================
general_ST_QCpipeline_R = function(counts, to_filter = FALSE)
{
  # ==========================================================================
  # STEP 1: load counts into a SingleCellExperiment
  # ==========================================================================
  sce = SingleCellExperiment(list(counts = counts))
  
  # ==========================================================================
  # STEP 2: compute per-cell QC metrics (before any filtering)
  # ==========================================================================
  message("############# QC b/ filtering #############")
  
  # recompute MT genes here, since some may have been filtered upstream
  mito.genes <- grep(pattern = "^MT-", x = rownames(counts(sce)), value = TRUE)
  
  sce <- addPerCellQC(sce, subsets = list(mito = mito.genes))
  
  # plot library size vs. number of genes detected, before any filtering
  plt_sum_detected = plotColData(sce, x = "sum", y = "detected")
  print(plt_sum_detected)
  
  #print(plotHighestExprs(sce, exprs_values = "counts", n = 50, feature_names_to_plot = "genes"))
  
  sce <- logNormCounts(sce)
  
  #vars <- getVarianceExplained(sce)
  #print(plotExplanatoryVariables(vars, ylim = c(0, 200)))
  
  # ==========================================================================
  # STEP 3: flag outlier cells (MAD-based thresholds) on library size, number of genes detected
  # and % mito and combine into a single pass/fail QC call
  # ==========================================================================
  
  libsize.drop_lower <- isOutlier(sce$total, nmads = 3, type = "lower", log = TRUE)
  libsize.drop_upper <- isOutlier(sce$total, nmads = 2, type = "higher", log = TRUE)
  
  feature.drop_lower <- isOutlier(sce$detected, nmads = 3, type = "lower", log = TRUE)
  feature.drop_upper <- isOutlier(sce$detected, nmads = 2, type = "higher", log = TRUE)
  
  mito.drop <- isOutlier(sce$subsets_mito_percent, nmads = 2.5, type = "higher")
  
  # a cell passes QC only if it's flagged as an outlier by none of the above
  pass <- !(libsize.drop_lower | libsize.drop_upper | feature.drop_lower | feature.drop_upper | mito.drop)
  sce$qc_pass <- factor(ifelse(test = pass,
                               yes = 'QC_pass',
                               no = 'QC_fail'))
  
  # ==========================================================================
  # STEP 4: diagnostic plots BEFORE filtering
  # (histograms with the MAD thresholds drawn as vertical lines)
  # ==========================================================================
  
  dt = colData(sce) %>% as.data.frame
  
  plt_hist_libsize_before_filtering = ggplot(dt, aes(x = total)) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.5) +
    xlab("Library size") +
    ylab("Number of cells") +
    geom_vline(xintercept = attr(libsize.drop_lower, 'thresholds')['lower'], col = "blue") +
    geom_vline(xintercept = attr(libsize.drop_upper, 'thresholds')['higher'], col = "blue") +
    scale_x_log10()
  print(plt_hist_libsize_before_filtering)
  
  plt_hist_Nexpressedgenes_before_filtering = ggplot(dt, aes(x = detected)) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.5) +
    xlab("Number of expressed genes") +
    ylab("Number of cells") +
    geom_vline(xintercept = attr(feature.drop_lower, 'thresholds')['lower'], col = "blue") +
    geom_vline(xintercept = attr(feature.drop_upper, 'thresholds')['higher'], col = "blue") +
    scale_x_log10()
  print(plt_hist_Nexpressedgenes_before_filtering)
  
  plt_hist_mitproportion_before_filtering = ggplot(dt, aes(x = subsets_mito_percent)) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.5) +
    xlab("Mitochondrial proportion (%)") +
    ylab("Number of cells") +
    geom_vline(xintercept = attr(mito.drop, 'thresholds')['higher'], col = "blue") +
    scale_x_log10()
  print(plt_hist_mitproportion_before_filtering)
  
  # library size vs. genes detected, colored by pass/fail - all cells, then
  # pass-only and fail-only separately, to visually sanity-check the QC call
  print(plotColData(sce, x = "sum", y = "detected", colour_by = 'qc_pass') +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position = "right", aspect.ratio = 1))
  
  print(plotColData(sce[, colData(sce)$qc_pass == "QC_pass"],
                    x = "sum", y = "detected", colour_by = 'qc_pass') +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position = "right", aspect.ratio = 1))
  
  print(plotColData(sce[, colData(sce)$qc_pass == "QC_fail"],
                    x = "sum", y = "detected", colour_by = 'qc_pass') +
          xlim(c(0, max(sce$sum))) +
          ylim(c(0, max(sce$detected))) +
          theme(legend.position = "right", aspect.ratio = 1))
  
  Sys.sleep(2) # let the previous plot finish rendering before the next message() below
  
  # ==========================================================================
  # STEP 5: apply filtering and diagnose the result
  # ==========================================================================
  message("############# QC a/ filtering #############")
  
  if (to_filter) {
    sce_filtered = sce[, pass]  # keep only cells that passed every QC check
  } else {
    sce_filtered = sce          # keep everything; QC metrics/plots are still informative
  }
  
  #print(plotHighestExprs(sce_filtered, exprs_values = "counts", n = 50, feature_names_to_plot = "genes"))
  
  # --- diagnostic plots AFTER filtering (same 3 histograms, no threshold lines needed) ---
  dt = colData(sce_filtered) %>% as.data.frame
  
  plt_hist_libsize_after_filtering = ggplot(dt, aes(x = total)) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.5) +
    xlab("Library size") +
    ylab("Number of cells") +
    scale_x_log10()
  print(plt_hist_libsize_after_filtering)
  
  plt_hist_Nexpressedgenes_after_filtering = ggplot(dt, aes(x = detected)) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.5) +
    xlab("Number of expressed genes") +
    ylab("Number of cells") +
    scale_x_log10()
  print(plt_hist_Nexpressedgenes_after_filtering)
  
  plt_hist_mitproportion_after_filtering = ggplot(dt, aes(x = subsets_mito_percent)) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.5) +
    xlab("Mitochondrial proportion (%)") +
    ylab("Number of cells") +
    scale_x_log10()
  print(plt_hist_mitproportion_after_filtering)
  
  Sys.sleep(2) # let the previous plot finish rendering before the next message() below
  
  # ==========================================================================
  # STEP 6: additional QC scatter/density plots (base R graphics, printed
  # directly to the active device rather than returned as ggplot objects)
  # ==========================================================================
  message("############# Additional QC metrics #############")
  
  # distribution of library sizes, with individual cells marked as a rug
  plot(density(sce_filtered$total / 1e3),
       xlab = "library size (thousands)",
       main = '')
  rug(sce_filtered$total / 1e3)
  
  # library size vs. mitochondrial proportion
  plot(y = sce_filtered$total / 1e3,
       x = sce_filtered$subsets_mito_percent,
       pch = 20,
       ylab = 'library size (thousands)',
       xlab = 'mitochondrial proportion (%)',
       ### col = ac('black', 0.5)
  )
  
  # library size vs. number of genes detected
  plot(y = sce_filtered$total / 1e3,
       x = sce_filtered$detected,
       pch = 20,
       ylab = 'library size (thousands)',
       xlab = 'number of genes',
       ### col = ac('black', 0.5)
  )
  
  # ==========================================================================
  # STEP 7: return the (filtered or unfiltered) SCE plus every diagnostic plot
  # ==========================================================================
  return(list(sce_filtered = sce_filtered,
              plt_sum_detected = plt_sum_detected,
              plt_hist_libsize_before_filtering = plt_hist_libsize_before_filtering,
              plt_hist_Nexpressedgenes_before_filtering = plt_hist_Nexpressedgenes_before_filtering,
              plt_hist_mitproportion_before_filtering = plt_hist_mitproportion_before_filtering,
              plt_hist_libsize_after_filtering = plt_hist_libsize_after_filtering,
              plt_hist_Nexpressedgenes_after_filtering = plt_hist_Nexpressedgenes_after_filtering,
              plt_hist_mitproportion_after_filtering = plt_hist_mitproportion_after_filtering))
}