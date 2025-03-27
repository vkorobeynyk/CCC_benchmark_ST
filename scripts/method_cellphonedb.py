from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import scanpy as sc
import pandas as pd

#############
### INPUT ###
#############
normalized_counts_path = snakemake.input["normalized_counts"]
microenvironment_path = snakemake.input["microenvironment"]
metadata_path = snakemake.input["metadata"]
cpdb_database_path = snakemake.input["cpdb_database"]
##############
### OUTPUT ###
##############
significant_interactions_path = snakemake.output["significant_interactions"]

#### Cellphonedb requires a file as input to the cpdb_statistical_analysis_method.
# I manually modified metadata file using the code below
#metadata= readRDS("data/processed/Visium_HPC_SJ/cellmetadata_Visium_HPC_SJ.RDS")$metadata
#colnames(metadata) = c("cell_type","barcode_sample")
#rownames(metadata) = metadata$barcode_sample
#write.table(metadata, "data/cpdbv5_extrafiles/metadata_for_cpdbv5.tsv" , sep = "\t")

#### Cellphonedb requires a file as input to the cpdb_statistical_analysis_method.
# I manually crated microenvironment file using the code below
#d = {"cell_type":["CT1","CT2"] , "microenvironment": ["env1","env1"]}
#microenvironment = pd.DataFrame(d)
#microenvironment.to_csv("data/Visium_HPC_SJ/microenvironment.tsv", sep = "\t", index = False)

out_path = 'results_cpdbv5_ToDelete/method1'

cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = cpdb_database_path,                 # mandatory: CellphoneDB database zip file.
  meta_file_path = metadata_path,                 # mandatory: tsv file defining barcodes to cell label.
  counts_file_path =  normalized_counts_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
  counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
  active_tfs_file_path = None,           # optional: defines cell types and their active TFs.
  microenvs_file_path = microenvironment_path,       # optional (default: None): defines cells per microenvironment.
  score_interactions = True,                       # optional: whether to score interactions or not. 
  iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
  threshold = 0,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
  threads = 5,                                     # number of threads to use in the analysis.
  debug_seed = 42,                                 # debug randome seed. To disable >=0.
  result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
  pvalue = 0.05,                                   # P-value threshold to employ for significance.
  subsampling = False,                             # To enable subsampling the data (geometri sketching).
  subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
  subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
  subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
  separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
  debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
  output_path = out_path,                          # Path to save results.
  output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
)

df = pd.DataFrame({"ligand_receptor" : cpdb_results["pvalues"]["interacting_pair"], "pval": cpdb_results["pvalues"]["CT1|CT2"]})
df = df[df["pval"] < 0.05]

df.to_csv(significant_interactions_path, sep = "\t")
