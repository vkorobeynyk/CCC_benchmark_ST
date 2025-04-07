import pandas as pd
import scanpy as sc
import liana as li
import numpy as np
import sys
import json
import yaml
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
plot_neighbors_path = snakemake.output["plot_neighbors"]

##############
### Params ###
##############
with open("config.yaml","r") as stream:
  config = yaml.safe_load(stream)

l_index = snakemake.params["l_index"]
dataset = snakemake.params["dataset"]

l = config["l_param"]["lianaP"][dataset][np.int64(l_index)]
#adata = sc.AnnData(pd.read_csv("output/STARmap_plus_HPC_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_4.tsv", sep="\t").T)
#cellmetadata_path = "output/STARmap_plus_HPC_semiSimulation_NB/simulated_cellmetadata_1_n_neigbors_4.json"

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# build the spatial graph with a selected bandwidth
li.ut.spatial_neighbors(adata, bandwidth=l, kernel='gaussian', set_diag=True)

#########################################
### Plot weights according to l param ###

df_neighbor_cells = pd.DataFrame(cellmetadata["neighbor_cells"])
Cell_OI = df_neighbor_cells.columns[2]
neighbor_cells = df_neighbor_cells.loc[:,Cell_OI]
# index of cell_OI 
cell_OI_index = np.where(cm["Cell_ID"] == Cell_OI)[0]

index_Positive_sp_connectivities_cell_OI = np.nonzero(adata.obsp["spatial_connectivities"][:,cell_OI_index] > 0.1)[0]
index_neighbor_cells = cm[cm["Cell_ID"].isin(list(neighbor_cells))].index

# Generate pplot
cm["spatial_connectivities"] = adata.obsp["spatial_connectivities"][:,cell_OI_index].A.flatten()

plt.scatter(cm["x"], cm["y"], 
            c=cm['spatial_connectivities'] , s = 5)
plt.colorbar()
plt.scatter(cm.loc[list(index_neighbor_cells),"x"], cm.loc[list(index_neighbor_cells),"y"], 
            c= "red", s = 5)
plt.scatter(cm.loc[cell_OI_index,"x"], cm.loc[cell_OI_index,"y"], 
            c= "green", s = 5)
            
plt.savefig(plot_neighbors_path, dpi = 200) 

# Bivariate Ligand-Receptor Relationships
lrdata = li.mt.bivariate(adata,
                      resource_name='cellchatdb', # NOTE: uses HUMAN gene symbols!
                      local_name='cosine', # Name of the function - currenty the other local metrics dont work/ dont change result at all
                      global_name="morans", # Name global function
                      n_perms=100, # Number of permutations to calculate a p-value
                      mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                      add_categories=True, # Whether to add local categories to the results
                      nz_prop=0, # Minimum expr. proportion for ligands/receptors and their subunits
                      seed=1,
                      use_raw=False,
                      verbose=True
)

#LRadata = out[1] # subset of adata to only LR only 
LRdata_df = lrdata.var.sort_values("morans_pvals", ascending=True) # extract df with interactions and statistics

# save data
LRdata_df = LRdata_df.loc[LRdata_df["morans_pvals"] < 0.05 , ["ligand","receptor","morans_pvals"]]
LRdata_df = LRdata_df.rename({"morans_pvals":"pval"},axis=1)
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]
LRdata_df.to_csv(significant_interactions_path, sep = "\t")
