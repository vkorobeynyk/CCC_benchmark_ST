import pandas as pd
import scanpy as sc
import spatialdm as sdm
import numpy as np
import math
import json
import matplotlib.pyplot as plt

#############
### INPUT ###
#############
processed_counts_path = snakemake.input["processed_counts"]
cellmetadata3_path = snakemake.input["cellmetadata3"]

##############
### OUTPUT ###
##############
cellmetadata4_path = snakemake.output["cellmetadata4"]
plot_neighbors_path = snakemake.output["plot_neighbors"]

adata = sc.AnnData(pd.read_csv(processed_counts_path, sep="\t").T)
with open(cellmetadata3_path,"r") as f:
  cellmetadata = json.load(f)

#processed_counts_path = "data/processed/STARmap_plus_HPC/processed_counts_STARmap_plus_HPC.tsv"
#cellmetadata2_path = "data/processed/STARmap_plus_HPC/cellmetadata2_STARmap_plus_HPC.json"
cm = pd.DataFrame(cellmetadata["metadata"])
df_neighbor_cells = pd.DataFrame(cellmetadata["neighbor_cells"])

# add x and y coordinates
adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# check that all elements in array are the same
def check(arr):
    return any(i != 6 for i in arr)
    
amount_neighbors_seen = [0]
l = 0.001
while check(amount_neighbors_seen):
    
    amount_neighbors_seen = []
    
    l = l+np.sqrt(l)
    # compute weights
    sdm.weight_matrix(adata, l=l, single_cell=False,n_neighbors = 300)
    
    for Cell_OI in df_neighbor_cells.columns:
        neighbor_cells = df_neighbor_cells[Cell_OI]
        # index of cell_OI 
        cell_OI_index = np.where(cm["Cell_ID"] == Cell_OI)[0]
        
        index_Positive_weights_cell_OI = np.nonzero(adata.obsp["weight"][:,cell_OI_index] > 0.0025)[0]
        index_neighbor_cells = cm[cm["Cell_ID"].isin(list(neighbor_cells))].index
            
        # if not all 6 possible neighbors are within the method sight, then break
        amount_neighbors_seen.append(len(set(index_neighbor_cells).intersection(set(index_Positive_weights_cell_OI))))
        
# transform all non-zero weight values to 1. 
# this is because the values of weights are so low (e-170) that the plot wont show them (and i didnt find other solution)
#adata.obsp['weight'][adata.obsp['weight'] > 0] = 1

# save figure showing sender , receiver and all cells with positive weights of last cell
# save figure showing sender , receiver and all cells with positive weights of last cell
plt.scatter(adata.obsm['spatial'][:,0], adata.obsm['spatial'][:,1], 
            c=adata.obsp['weight'].A[cell_OI_index] , s = 5)
plt.colorbar()
plt.scatter(adata.obsm['spatial'][index_neighbor_cells,0], adata.obsm['spatial'][index_neighbor_cells,1], 
            c= "red", s = 5)
plt.scatter(adata.obsm['spatial'][cell_OI_index,0], adata.obsm['spatial'][cell_OI_index,1], 
            c= "green", s = 5)
            
plt.savefig(plot_neighbors_path, dpi = 200) 

cellmetadata["spatialWeight_spatialdm"] = l

# save metadata as json file
with open(cellmetadata4_path, "w") as outfile: 
  json.dump(cellmetadata, outfile)
