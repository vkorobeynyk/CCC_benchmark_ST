import pandas as pd
import scanpy as sc
import liana as li
import numpy as np
import sys
import json
import matplotlib.pyplot as plt

#############
### INPUT ###
#############
processed_counts_path = snakemake.input["processed_counts"]
cellmetadata2_path = snakemake.input["cellmetadata2"]

##############
### OUTPUT ###
##############
cellmetadata3_path = snakemake.output["cellmetadata3"]
plot_neighbors_path = snakemake.output["plot_neighbors"]

adata = sc.AnnData(pd.read_csv(processed_counts_path, sep="\t").T)
with open(cellmetadata2_path,"r") as f:
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
l = 0.1
while check(amount_neighbors_seen):
    
    amount_neighbors_seen = []
    
    l = l + 0.1
    # compute weights
    li.ut.spatial_neighbors(adata, bandwidth=l, kernel='gaussian', set_diag=True,max_neighbours=300)
    
    for Cell_OI in df_neighbor_cells.columns:
        neighbor_cells = df_neighbor_cells[Cell_OI]
        # index of cell_OI 
        cell_OI_index = np.where(cm["Cell_ID"] == Cell_OI)[0]
        
        index_Positive_sp_connectivities_cell_OI = np.nonzero(adata.obsp["spatial_connectivities"][:,cell_OI_index] > 0.1)[0]
        index_neighbor_cells = cm[cm["Cell_ID"].isin(list(neighbor_cells))].index
            
        # if not all 6 possible neighbors are within the method sight, then break
        amount_neighbors_seen.append(len(set(index_neighbor_cells).intersection(set(index_Positive_sp_connectivities_cell_OI))))
        

# save figure showing sender , receiver and all cells with positive weights of last cell
cm["spatial_connectivities"] = adata.obsp["spatial_connectivities"][:,cell_OI_index].A.flatten()

plt.scatter(cm["x"], cm["y"], 
            c=cm['spatial_connectivities'] , s = 5)
plt.colorbar()
plt.scatter(cm.loc[list(index_neighbor_cells),"x"], cm.loc[list(index_neighbor_cells),"y"], 
            c= "red", s = 5)
plt.scatter(cm.loc[cell_OI_index,"x"], cm.loc[cell_OI_index,"y"], 
            c= "green", s = 5)
            
plt.savefig(plot_neighbors_path, dpi = 200) 

cellmetadata["spatialWeight_lianaPlus"] = l

# save metadata as json file
with open(cellmetadata3_path, "w") as outfile: 
  json.dump(cellmetadata, outfile)
