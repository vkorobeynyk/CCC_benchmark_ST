import pandas as pd
import scanpy as sc
import spatialdm as sdm
import numpy as np
import math
import json
import matplotlib.pyplot as plt
import yaml

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
plot_neighbors_path = snakemake.output["plot_neighbors"]

##############
### Params ###
##############
with open("config.yaml","r") as stream:
  config = yaml.safe_load(stream)

l_index = snakemake.params["l_index"]
dataset = snakemake.params["dataset"]
print(config)
l = config["l_param"]["spatialdm"][dataset][np.int64(l_index)]

#adata = sc.AnnData(pd.read_csv("output/STARmap_plus_HPC_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_4.tsv", sep="\t").T)
#cellmetadata_path = "output/STARmap_plus_HPC_semiSimulation_NB/simulated_cellmetadata_1_n_neigbors_4.json"

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# create weight matrix by rbf kernel
sdm.weight_matrix(adata, l=l, single_cell=False)

#########################################
### Plot weights according to l param ###

# transform all non-zero weight values to 1. 
# this is because the values of weights are so low (e-170) that the plot wont show them (and i didnt find other solution)
adata.obsp['weight'][adata.obsp['weight'] > 0] = 1

df_neighbor_cells = pd.DataFrame(cellmetadata["neighbor_cells"])
Cell_OI = df_neighbor_cells.columns[2]
neighbor_cells = df_neighbor_cells.loc[:,Cell_OI]
# index of cell_OI 
cell_OI_index = np.where(cm["Cell_ID"] == Cell_OI)[0]
index_neighbor_cells = cm[cm["Cell_ID"].isin(list(neighbor_cells))].index

# save figure showing sender , receiver and all cells with positive weights of last cell
plt.scatter(adata.obsm['spatial'][:,0], adata.obsm['spatial'][:,1], 
            c=adata.obsp['weight'].A[cell_OI_index] , s = 5)
plt.colorbar()
plt.scatter(adata.obsm['spatial'][index_neighbor_cells,0], adata.obsm['spatial'][index_neighbor_cells,1], 
            c= "red", s = 5)
plt.scatter(adata.obsm['spatial'][cell_OI_index,0], adata.obsm['spatial'][cell_OI_index,1], 
            c= "black", s = 5)
plt.savefig(plot_neighbors_path, dpi = 200)

# extract LR
sdm.extract_lr(adata, 'human', min_cell=0) # uses cellchatdb by default

# global Moran selection
sdm.spatialdm_global(adata, 1000, specified_ind=None, method='both', nproc=1)

# select significant pairs
sdm.sig_pairs(adata, method='permutation', fdr=True, threshold=0.1)     

LRdata_df = adata.uns['global_res'].sort_values("perm_pval", ascending=True)
LRdata_df = LRdata_df.loc[:  , ["Ligand0","Receptor0","perm_pval"]]
LRdata_df = LRdata_df.rename({"Ligand0":"ligand" , "Receptor0":"receptor", "perm_pval":"statistics"},axis=1)
LRdata_df = LRdata_df.sort_values("statistics") # sort importance column on ascending order
LRdata_df["significant"] = LRdata_df["statistics"] < 0.05
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]

LRdata_df.to_csv(significant_interactions_path, sep = "\t")
