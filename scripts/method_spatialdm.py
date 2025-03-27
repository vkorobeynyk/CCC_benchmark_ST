import pandas as pd
import scanpy as sc
import spatialdm as sdm
import numpy as np
import math
import json
from rds2py import read_rds
import matplotlib.pyplot as plt

#############
### INPUT ###
#############
normalized_counts_path = snakemake.input["normalized_counts"]
cellmetadata_path = snakemake.input["cellmetadata_post_simulation"]

adata = sc.AnnData(pd.read_csv(normalized_counts_path, sep="\t").T)
with open(cellmetadata_path,"r") as f:
  cellmetadata = json.load(f)
  
cm = pd.DataFrame(cellmetadata["metadata"])

##############
### OUTPUT ###
##############
significant_interactions_path = snakemake.output["significant_interactions"]

##############
### Params ###
##############
l = cellmetadata["spatialWeight_spatialdm"]

#adata = sc.AnnData(pd.read_csv("output/STARmap_plus_HPC_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_4.tsv", sep="\t").T)
#cellmetadata_path = "output/STARmap_plus_HPC_semiSimulation_NB/simulated_cellmetadata_1_n_neigbors_4.json"

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# create weight matrix by rbf kernel
sdm.weight_matrix(adata, l=l[0], single_cell=False)

sdm.extract_lr(adata, 'human', min_cell=0) # uses cellchatdb by default

# global Moran selection
sdm.spatialdm_global(adata, 1000, specified_ind=None, method='both', nproc=1)

# select significant pairs
sdm.sig_pairs(adata, method='permutation', fdr=True, threshold=0.1)     

LRdata_df = adata.uns['global_res'].sort_values("perm_pval", ascending=True)
LRdata_df = LRdata_df.loc[LRdata_df["perm_pval"] < 0.05 , ["Ligand0","Receptor0","perm_pval"]]
LRdata_df = LRdata_df.rename({"Ligand0":"ligand" , "Receptor0":"receptor", "perm_pval":"pval"},axis=1)
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]

LRdata_df.to_csv(significant_interactions_path, sep = "\t")
