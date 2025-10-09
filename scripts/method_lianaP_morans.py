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
LR_database_path = snakemake.params["LR_database"]

l = config["l_param"]["lianaP"][dataset][np.int64(l_index)]
#adata = sc.AnnData(pd.read_csv("output/MERFISH_mColon_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_2_indexLR_1.tsv", sep="\t").T)
#cellmetadata_path = "output/STARmap_plus_HPC_semiSimulation_NB/simulated_cellmetadata_1_n_neigbors_4.json"

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# build the spatial graph with a selected bandwidth
li.ut.spatial_neighbors(adata, bandwidth=l, kernel='gaussian', set_diag=True)

#########################################
### Plot weights according to l param ###

all_sender_cells = cm.loc[cm["Celltype_updated"] == "CT1_signalAdded",]["Cell_ID"].tolist()
index_sender_Cells = cm.loc[cm["Celltype_updated"] == "CT1_signalAdded",]["Cell_ID"].index.tolist()
all_receiver_cells = cm.loc[cm["Celltype_updated"] == "CT2_signalAdded",]["Cell_ID"].tolist()
index_receiver_cells = cm.loc[cm["Celltype_updated"] == "CT2_signalAdded",]["Cell_ID"].index.tolist()

# get spatial connectivities
cm["spatial_connectivities"] = pd.DataFrame.sparse.from_spmatrix(adata.obsp["spatial_connectivities"]).loc[index_sender_Cells,].max(axis=0)
            
plt.scatter(cm["x"], cm["y"], 
            c=cm['spatial_connectivities'] , s = 5)
plt.colorbar()
plt.scatter(cm.loc[cm["Cell_ID"].isin(all_sender_cells),"x"], cm.loc[cm["Cell_ID"].isin(all_sender_cells),"y"], 
            c= "#990099", s = 5)
plt.scatter(cm.loc[cm["Cell_ID"].isin(all_receiver_cells),"x"], cm.loc[cm["Cell_ID"].isin(all_receiver_cells),"y"], 
            c= "#0000FF", s = 5)
plt.xlabel("x_coord_um")
plt.ylabel("y_coord_um")
plt.title("purple - sender cells | blue - receiver cells | orange - cells seen by method | color bar - spatial connectivity values",
          fontsize = 2)
            
plt.savefig(plot_neighbors_path, dpi = 200) 

# get LR database
LR_database = pd.read_csv(LR_database_path, sep=" ")

# Bivariate Ligand-Receptor Relationships
liana = li.mt.bivariate(adata,
                      resource=LR_database,
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

# extract info and save data
LRdata_df = liana.var
LRdata_df = LRdata_df.loc[:  , ["ligand","receptor","morans_pvals"]]
LRdata_df["significant"] = LRdata_df["morans_pvals"] < 0.05
LRdata_df = LRdata_df.rename({"morans_pvals":"statistics"},axis=1)
LRdata_df = LRdata_df.sort_values("statistics") # sort importance column on ascending order
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]
LRdata_df.to_csv(significant_interactions_path, sep = "\t")
