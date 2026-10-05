# ============================================================================
# Baseline method: uses a simple Seurat Wilcoxon test to see how robustly a
# generic differential-expression approach detects the simulated L-R
# interaction.
# ============================================================================

suppressMessages({
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(jsonlite)
  library(magrittr)
  library(purrr)
  source("scripts/helper_functions.R")
})

# Fail early with a clear error if a required input/output is missing
if (is.null(snakemake@input[["normalized_counts"]]) | is.null(snakemake@input[["cellmetadata_post_simulation"]])
    | is.null(snakemake@output[["significant_interactions"]])) {
  stop("Argument_name needs to be specified, but is missing.n", call. = FALSE)
}

# ==============================================================================
# STEP 0: Snakemake I/O and params
# ==============================================================================

# input files
normalized_counts_path = snakemake@input[["normalized_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata_post_simulation"]]

# output files
significant_interactions_path = snakemake@output[["significant_interactions"]]

# params
dataset = snakemake@params["dataset"] %>% as.character

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)

# ==============================================================================
# STEP 1: load data
# ==============================================================================

inflated_counts = read.csv(normalized_counts_path, sep = "\t") %>% as.matrix
cellmetadata = read_json(path = cellmetadata_path)

# transform the JSON list into individual dataframes
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

# Create Seurat object (will contain only normalized counts)
SO_obj = CreateSeuratObject(counts = inflated_counts, meta.data = cellmetadata$metadata)
SO_obj@assays$RNA$data = inflated_counts

# ==============================================================================
# STEP 2: expand the LR database to single ligand-receptor gene pairs
#
# Seurat only does pairwise comparisons without subunit info, so cases like
# L-R1-R2 aren't directly possible. To work around this, we split them into
# L-R1 and L-R2 and average the results afterwards.
# ==============================================================================

# Filter LR_database and the Seurat object down to their common ligand/receptor genes.
# Transform LR database from L-R1-R2 to L-R1 | L-R2 rows.
lst = split(LR_database, seq_len(nrow(LR_database)))

LR_database = lapply(lst, function(x) {
  grid = x[c("ligand", "receptor")] %>%
    str_split(., "_") %>%
    expand.grid() %>%
    mutate(ligand_receptor = x$ligand_receptor)
  
  colnames(grid)[1:2] = c("ligand", "receptor")
  all_present = (grid %>% select("ligand", "receptor") %>% unlist) %in% rownames(SO_obj) %>%
    all
  
  if (all_present) {return(grid)}
}) %>% Reduce(rbind.data.frame, .)

# remove duplicated entries: some expanded interactions coincide (e.g. WNT2B_FZD7_LRP6 and WNT2B_FZD4_LRP6)
n = LR_database %>%
  select(ligand, receptor) %>%
  apply(., 1, function(x) {str_c(x, collapse = "_")}) %>%
  duplicated %>%
  which

LR_database = LR_database[-n, ]

# filter the Seurat object down to only the L/R genes
SO_obj = SO_obj[c(LR_database$ligand, LR_database$receptor) %>% unique %>% as.character, ]

# ==============================================================================
# STEP 3: Wilcoxon test - find ligands over-expressed in Sender and
# receptors over-expressed in Receiver
#
# Some genes may be dropped from the output due to min.cells.feature/min.cells.group.
# ==============================================================================

Idents(SO_obj) = SO_obj$Celltype
markers = FindAllMarkers(SO_obj, slot = "data", test.use = "wilcox")

# Initialize as empty dataframes with the correct column names if no markers are found
if (nrow(markers) == 0) {
  markers_ligands = markers
  markers_receptors = markers
} else {
  markers_ligands = markers %>% filter(cluster == "Sender") %>% filter(avg_log2FC > 0)
  markers_receptors = markers %>% filter(cluster == "Receiver") %>% filter(avg_log2FC > 0)
}

# ==============================================================================
# STEP 4: recombine ligand + receptor(subunit) marker results per LR pair
# ==============================================================================

seurat_LRR_averaged_out = list()

if (nrow(markers_ligands) > 0 & nrow(markers_receptors) > 0) {
  
  for (entry in unique(LR_database$ligand_receptor)) {
    name = unlist(str_split(entry, "_"))
    ligand = name[1]
    
    # handle cases with 1 or 2 receptor subunits (L_R or L_R1_R2)
    if (length(name) %in% c(2, 3)) {
      receptors = name[2:length(name)]
      
      # subset markers for this specific pair
      subset_L = markers_ligands %>% filter(gene == ligand)
      subset_R = markers_receptors %>% filter(gene %in% receptors)
      
      # only keep it if ALL components were found (1 ligand, and as many
      # receptor subunits as the database string lists)
      if (nrow(subset_L) == 1 && nrow(subset_R) == (length(name) - 1)) {
        
        tmp_markers = rbind(subset_L, subset_R)
        
        averaged_results = tmp_markers %>%
          dplyr::summarise(
            across(where(is.character), ~ data.table::first(.)),
            across(where(is.factor), ~ data.table::first(.)),
            across(where(is.numeric), ~ mean(., na.rm = TRUE))
          )
        
        seurat_LRR_averaged_out[[entry]] = averaged_results
      }
    } else {
      message(paste("Skipping complex pair:", entry))
    }
  }
}

if (length(seurat_LRR_averaged_out) > 0) {
  seurat_LRR_averaged_out = do.call(rbind.data.frame, seurat_LRR_averaged_out) %>%
    mutate(
      ligand_receptor = rownames(.),
      significant = p_val_adj < 0.05,
      statistics = p_val_adj
    ) %>%
    select(p_val, p_val_adj, ligand_receptor, significant, statistics) %>%
    arrange(statistics)
} else {
  # default output when no LR pairs were significantly co-expressed
  seurat_LRR_averaged_out = data.frame(
    p_val = 1,
    p_val_adj = 1,
    ligand_receptor = "none_detected",
    significant = FALSE,
    statistics = 1
  )
}

# metadata columns (this method doesn't use a spatial radius, so these are
# fixed values rather than computed values, unlike the spatial methods)
seurat_LRR_averaged_out$ratio_Receiver_seen_byMethod = 100
seurat_LRR_averaged_out$average_cells_perSender_seen_byMethod = FALSE

# ==============================================================================
# STEP 5: save results
# ==============================================================================

write.table(seurat_LRR_averaged_out, significant_interactions_path, sep = "\t")