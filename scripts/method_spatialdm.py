import pandas as pd
import scanpy as sc
import spatialdm as sdm
import numpy as np
import math
import json
import matplotlib.pyplot as plt
import yaml
import inspect
import re

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
l = config["l_param"]["spatialdm"][dataset][np.int64(l_index)]

#adata = sc.AnnData(pd.read_csv("output/STARmap_plus_HPC_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_4.tsv", sep="\t").T)
#cellmetadata_path = "output/STARmap_plus_HPC_semiSimulation_NB/simulated_cellmetadata_1_n_neigbors_4.json"

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# spatiualDM download the interaction files from figshare.
# there is an issue in source code with parsing, so I download the files manually and change source code of the function here.
# Get the original source
source = inspect.getsource(sdm.extract_lr)

# Replace the broken URL line with your local path logic
# (You'll need to check the exact variable name in the source)
new_source = re.sub("https://figshare.com/ndownloader/files/36638940", "data/complex_input_CellChatDB.csv", source)
new_source = re.sub("https://figshare.com/ndownloader/files/36638943", "data/interaction_input_CellChatDB.csv", new_source)

# Inject it into the module dictionary
exec(new_source, sdm.__dict__)

# create weight matrix by rbf kernel
sdm.weight_matrix(adata, l=l, single_cell=False)

#########################################
### Plot weights according to l param ###

# transform all non-zero weight values to 1. 
# this is because the values of weights are so low (e-170) that the plot wont show them (and i didnt find other solution)
#adata.obsp['weight'][adata.obsp['weight'] > 0] = 1

df_neighbor_cells = pd.DataFrame(cellmetadata["neighbor_cells"])
Cell_OI = df_neighbor_cells.columns[2]
neighbor_cells = df_neighbor_cells.loc[:,Cell_OI]
# index of cell_OI 
cell_OI_index = np.where(cm["Cell_ID"] == Cell_OI)[0]
index_neighbor_cells = cm[cm["Cell_ID"].isin(list(neighbor_cells))].index

# save figure showing sender , receiver and all cells with positive weights of last cell
plt.scatter(adata.obsm['spatial'][:,0], adata.obsm['spatial'][:,1], 
            c=adata.obsp['weight'].toarray()[cell_OI_index] , s = 5)
plt.colorbar()
plt.scatter(adata.obsm['spatial'][index_neighbor_cells,0], adata.obsm['spatial'][index_neighbor_cells,1], 
            c= "red", s = 5)
plt.scatter(adata.obsm['spatial'][cell_OI_index,0], adata.obsm['spatial'][cell_OI_index,1], 
            c= "black", s = 5)
plt.savefig(plot_neighbors_path, dpi = 200)

# extract LR
sdm.extract_lr(adata, 'human', min_cell=0) # uses cellchatdb by default

# global Moran selection
sdm.spatialdm_global(adata, n_perm=250, specified_ind=None, method='both', nproc=1)

# select significant pairs
sdm.sig_pairs(adata, method='permutation', fdr=True, threshold=0.1)     

all_sender_cells = cm.loc[cm["Celltype_updated"] == "CT1_signalAdded",]["Cell_ID"].tolist()
index_sender_Cells = cm.loc[cm["Celltype_updated"] == "CT1_signalAdded",]["Cell_ID"].index.tolist()
all_receiver_cells = cm.loc[cm["Celltype_updated"] == "CT2_signalAdded",]["Cell_ID"].tolist()
index_receiver_cells = cm.loc[cm["Celltype_updated"] == "CT2_signalAdded",]["Cell_ID"].index.tolist()

# get spatial weights
cm["weight"] = pd.DataFrame.sparse.from_spmatrix(adata.obsp["weight"]).loc[index_sender_Cells,].max(axis=0)

# how many CT2 cells are seen by method
amount_CT2_seen_byMethod = (cm[cm["Cell_ID"].isin(all_receiver_cells)]["weight"]!= 0).sum() 

# FCsender > FCreceiver -> how many receivers are seen by method?
# FCsender < FCreceiver -> do all senders see 1 receiver?
# FCsender == FCreceiver -> are all receiver seen by CT1 and method
if amount_CT2_seen_byMethod != 0:
    if len(all_sender_cells) > len(all_receiver_cells):
        ratio_CT2_seen_byMethod = amount_CT2_seen_byMethod / len(all_receiver_cells) * 100
    elif len(all_sender_cells) < len(all_receiver_cells):
        ratio_CT2_seen_byMethod = len(all_sender_cells) / amount_CT2_seen_byMethod * 100
        # if there are more than 1 receiver per sender
        #if amount_CT2_seen_byMethod > len(all_sender_cells):
        #    ratio_CT2_seen_byMethod = 100
    elif len(all_sender_cells) == len(all_receiver_cells):
        ratio_CT2_seen_byMethod = amount_CT2_seen_byMethod / len(all_receiver_cells) * 100
else:
    ratio_CT2_seen_byMethod = 0
    
# average cells that each CT1 has that are seen by method
average_cells_perCT1_seen_byMethod = (cm["weight"]!=0).sum() / len(all_sender_cells)


LRdata_df = adata.uns['global_res'].sort_values("perm_pval", ascending=True)
LRdata_df = LRdata_df.loc[:  , ["Ligand0","Receptor0","Receptor1","perm_pval"]]
LRdata_df = LRdata_df.rename({"Ligand0":"ligand" , "perm_pval":"statistics"},axis=1)
LRdata_df = LRdata_df.sort_values("statistics") # sort importance column on ascending order
LRdata_df["significant"] = LRdata_df["statistics"] < 0.05
LRdata_df["receptor"] = LRdata_df["Receptor0"] + "_" + LRdata_df["Receptor1"].fillna("")
LRdata_df["receptor"] = LRdata_df["receptor"].str.rstrip("_")
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]
LRdata_df["ratio_CT2_seen_byMethod"] = ratio_CT2_seen_byMethod
LRdata_df["average_cells_perCT1_seen_byMethod"] = average_cells_perCT1_seen_byMethod

LRdata_df.to_csv(significant_interactions_path, sep = "\t")
