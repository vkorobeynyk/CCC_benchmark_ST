import pandas as pd
import scanpy as sc
from liana.method import lrMistyData
from liana.method.sp import LinearModel
import numpy as np
import sys
import json
import yaml
import matplotlib.pyplot as plt

# Docs -> https://liana-py.readthedocs.io/en/latest/notebooks/misty.html#Ligand-Receptor-Misty

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

# get LR database
LR_database = pd.read_csv(LR_database_path, sep=" ")

# Build LR Misty object:
misty = lrMistyData(adata,
                      resource=LR_database,
                      nz_threshold=0, # this make paraview also see direct neighboorhods (which is supposed to be seen by only justaview)
                      cutoff=0.01, # doesnt influence output too much
                      spatial_key='spatial',
                      kernel="gaussian",
                      bandwidth=l,
                      use_raw=False,
                      verbose=True
)

misty(bypass_intra=True, model=LinearModel, verbose=True)

#########################################
### Plot weights according to l param ###


all_sender_cells = cm.loc[cm["Celltype_updated"] == "CT1_signalAdded",]["Cell_ID"].tolist()
index_sender_Cells = cm.loc[cm["Celltype_updated"] == "CT1_signalAdded",]["Cell_ID"].index.tolist()
all_receiver_cells = cm.loc[cm["Celltype_updated"] == "CT2_signalAdded",]["Cell_ID"].tolist()
index_receiver_cells = cm.loc[cm["Celltype_updated"] == "CT2_signalAdded",]["Cell_ID"].index.tolist()

# get spatial connectivities compputd by misty for all sender cells
cm["spatial_connectivities"] = pd.DataFrame.sparse.from_spmatrix(misty["extra"].obsp["spatial_connectivities"]).loc[index_sender_Cells,].max(axis=0)
            
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



# save data
LRdata_df = misty.uns["interactions"]
# importances relate to positive/negative correlation
# in cases where we have a disbalanced amount of ligand/receptor, we have negative correlation, that is still an interesting observation to take into account
# so we convert all importances to positive for ranking purposes
LRdata_df["importances"] = LRdata_df["importances"].abs() 
LRdata_df["statistics"] = LRdata_df["importances"] # for ranking inflated LR we are only interested in communication between sender and environment and not environment to sender
LRdata_df = LRdata_df.sort_values("statistics", ascending=False) # sort importance column on descending order (higher importance first)
LRdata_df["ligand_receptor"] = LRdata_df["predictor"] + "_" + LRdata_df["target"] # intra view we have receptors and para view ligands. Thats why its predictor_target
LRdata_df["significant"] = LRdata_df["importances"] > 2 # importance of 2 was used in the paper as a filtering criteria

LRdata_df.to_csv(significant_interactions_path, sep = "\t")
