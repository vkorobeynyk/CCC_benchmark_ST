"""
method_cellphonedb.py

Runs CellPhoneDB (v5, statistical analysis method) on the semi-simulated
spatial dataset and extracts the significance of the simulated
Sender -> Receiver ligand-receptor interaction.
"""

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import scanpy as sc
import pandas as pd

# ==============================================================================
# STEP 0: Snakemake I/O
# ==============================================================================

# input files
normalized_counts_path = snakemake.input["normalized_counts"]
microenvironment_path = snakemake.input["microenvironment"]
cellmetadata_path = snakemake.input["cellmetadata_post_simulation"]
cpdb_database_path = snakemake.input["cpdb_database"]

# output files
significant_interactions_path = snakemake.output["significant_interactions"]

# ==============================================================================
# STEP 1: run CellPhoneDB
#
# CellPhoneDB requires a metadata file as input to cpdb_statistical_analysis_method.
# That file is generated in processing_semisimulation_nbsampling.R.
# The microenvironment file was created manually once, using:
#   d = {"cell_type": ["Sender", "Receiver"], "microenvironment": ["env1", "env1"]}
#   microenvironment = pd.DataFrame(d)
#   microenvironment.to_csv("data/Visium_HPC_SJ/microenvironment.tsv", sep="\t", index=False)
# ==============================================================================

out_path = 'results_cpdbv5_ToDelete/method1'

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path=cpdb_database_path,           # mandatory: CellPhoneDB database zip file.
    meta_file_path=cellmetadata_path,            # mandatory: tsv file defining barcodes to cell label.
    counts_file_path=normalized_counts_path,     # mandatory: normalized count matrix - path or in-memory AnnData object.
    counts_data='hgnc_symbol',                   # defines the gene annotation used in the count matrix.
    active_tfs_file_path=None,                   # optional: defines cell types and their active TFs.
    microenvs_file_path=microenvironment_path,   # optional (default: None): defines cells per microenvironment.
    score_interactions=True,                     # optional: whether to score interactions.
    iterations=1000,                             # number of shufflings performed in the analysis.
    threshold=0,                                 # min % of cells expressing a gene for it to be used in the analysis.
    threads=1,                                   # number of threads to use.
    debug_seed=42,                               # debug random seed (disable with a value >= 0... actually: to disable, use -1).
    result_precision=3,                          # rounding for the mean values in significant_means.
    pvalue=1,                                    # p-value threshold used for significance.
    subsampling=False,                           # whether to subsample the data (geometric sketching).
    subsampling_log=False,                       # (mandatory) enable subsampling log1p for non-log-transformed inputs.
    subsampling_num_pc=100,                      # number of components to subsample via geometric sketching (default: 100).
    subsampling_num_cells=1000,                  # number of cells to subsample (default: 1/3 of the dataset).
    separator='|',                               # string used to separate cells in the results, e.g. "cellA|cellB".
    debug=False,                                 # save all intermediate tables used during analysis (pkl format).
    output_path=out_path,                        # path to save results.
    output_suffix=None                           # replaces the timestamp in output filenames with a custom string (default: None).
)

# ==============================================================================
# STEP 2: clean up multi-subunit interaction names
#
# CellPhoneDB reports multi-subunit interactions as, e.g., "FN1_integrin_aVb1_complex".
# Replace those with the explicit gene list, e.g. "FN1_ITGB1_ITGAV".
# ==============================================================================

cpdb_results["pvalues"] = cpdb_results["pvalues"].reset_index(drop=True)
for row, pair in cpdb_results["pvalues"]["interacting_pair"].items():
    if "complex" in str(pair):
        interaction_id = cpdb_results["pvalues"]["id_cp_interaction"].iloc[row]
        interaction_id_all_genes = cpdb_results["deconvoluted"][cpdb_results["deconvoluted"]["id_cp_interaction"] == interaction_id]
        interaction_id_joined = "_".join(interaction_id_all_genes["gene_name"].astype(str).tolist())

        cpdb_results["pvalues"]["interacting_pair"][row] = interaction_id_joined

# ==============================================================================
# STEP 3: extract the Sender -> Receiver p-values and save results
# ==============================================================================

LRdata_df = pd.DataFrame({"ligand_receptor": cpdb_results["pvalues"]["interacting_pair"], "pval": cpdb_results["pvalues"]["Sender|Receiver"]})
LRdata_df["significant"] = LRdata_df["pval"] < 0.05
LRdata_df = LRdata_df.rename({"pval": "statistics"}, axis=1)
LRdata_df = LRdata_df.sort_values("statistics")  # sort importance column in ascending order

# metadata columns (this method doesn't use a spatial radius, so these are
# fixed placeholders rather than computed values, unlike the spatial methods)
LRdata_df["ratio_Receiver_seen_byMethod"] = 100
LRdata_df["average_cells_perSender_seen_byMethod"] = False

LRdata_df.to_csv(significant_interactions_path, sep="\t")
