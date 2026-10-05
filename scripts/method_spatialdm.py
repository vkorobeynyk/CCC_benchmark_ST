"""
method_spatialdm.py

Runs SpatialDM on the semi-simulated spatial dataset and extracts the
significance of the simulated Sender -> Receiver ligand-receptor interaction.
"""

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

# ==============================================================================
# STEP 0: Snakemake I/O and params
# ==============================================================================

# input files
normalized_counts_path = snakemake.input["normalized_counts"]
cellmetadata_path = snakemake.input["cellmetadata_post_simulation"]

adata = sc.AnnData(pd.read_csv(normalized_counts_path, sep="\t").T)
with open(cellmetadata_path, "r") as f:
  cellmetadata = json.load(f)

cm = pd.DataFrame(cellmetadata["metadata"])

# output files
significant_interactions_path = snakemake.output["significant_interactions"]
plot_neighbors_path = snakemake.output["plot_neighbors"]

# params
with open("config.yaml", "r") as stream:
  config = yaml.safe_load(stream)

radius_index = snakemake.params["radius_index"]
dataset = snakemake.params["dataset"]
l = config["radius_param"]["spatialdm"][dataset][np.int64(radius_index)]

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# ==============================================================================
# STEP 1: build the spatial weight matrix (RBF kernel)
#
# SpatialDM normally downloads its interaction files from figshare. There's a
# parsing issue in the upstream source for that, so we download the files
# manually and patch the function's source at runtime to point to our local
# copies.
# ==============================================================================

source = inspect.getsource(sdm.extract_lr)

# Replace the broken URLs with local file paths
new_source = re.sub("https://figshare.com/ndownloader/files/36638940", "data/SPATIALDM_REQ_complex_input_CellChatDB.csv", source)
new_source = re.sub("https://figshare.com/ndownloader/files/36638943", "data/SPATIALDM_REQ_interaction_input_CellChatDB.csv", new_source)

# inject the patched function back into the module
exec(new_source, sdm.__dict__)

sdm.weight_matrix(adata, l=l, single_cell=False)

# ==============================================================================
# STEP 2: diagnostic plot of a single Sender cell's spatial weights
# ==============================================================================

# weight values can be so small (e.g. ~1e-170) that
# they don't show up in the plot. Left here as a possible fix if needed.
# adata.obsp['weight'][adata.obsp['weight'] > 0] = 1

df_neighbor_cells = pd.DataFrame(cellmetadata["neighbor_cells"])
Cell_OI = df_neighbor_cells.columns[2]
neighbor_cells = df_neighbor_cells.loc[:, Cell_OI]
# index of Cell_OI
cell_OI_index = np.where(cm["Cell_ID"] == Cell_OI)[0]
index_neighbor_cells = cm[cm["Cell_ID"].isin(list(neighbor_cells))].index

# save a figure showing Sender, Receiver, and all cells with positive weight relative to the last cell
plt.scatter(adata.obsm['spatial'][:, 0], adata.obsm['spatial'][:, 1],
            c=adata.obsp['weight'].toarray()[cell_OI_index], s=5)
plt.colorbar()
plt.scatter(adata.obsm['spatial'][index_neighbor_cells, 0], adata.obsm['spatial'][index_neighbor_cells, 1],
            c="red", s=5)
plt.scatter(adata.obsm['spatial'][cell_OI_index, 0], adata.obsm['spatial'][cell_OI_index, 1],
            c="black", s=5)
plt.savefig(plot_neighbors_path, dpi=200)

# ==============================================================================
# STEP 3: run SpatialDM (LR extraction, global Moran selection, significance test)
# ==============================================================================

sdm.extract_lr(adata, 'human', min_cell=0)  # extract LR pairs (uses CellChatDB by default)
sdm.spatialdm_global(adata, n_perm=250, specified_ind=None, method='both', nproc=1)
sdm.sig_pairs(adata, method='permutation', fdr=True, threshold=0.1)

# ==============================================================================
# STEP 4: neighborhood-coverage metrics
# ==============================================================================

all_sender_cells = cm.loc[cm["Celltype_updated"] == "Sender_signalAdded", ]["Cell_ID"].tolist()
index_sender_Cells = cm.loc[cm["Celltype_updated"] == "Sender_signalAdded", ]["Cell_ID"].index.tolist()
all_receiver_cells = cm.loc[cm["Celltype_updated"] == "Receiver_signalAdded", ]["Cell_ID"].tolist()
index_receiver_cells = cm.loc[cm["Celltype_updated"] == "Receiver_signalAdded", ]["Cell_ID"].index.tolist()

# get spatial weights
cm["weight"] = pd.DataFrame.sparse.from_spmatrix(adata.obsp["weight"]).loc[index_sender_Cells, ].max(axis=0)

# how many Receiver cells are seen by the method
amount_Receiver_seen_byMethod = (cm[cm["Cell_ID"].isin(all_receiver_cells)]["weight"] != 0).sum()

# PCE_Sender > PCE_Receiver -> how many Receivers are seen by the method?
# PCE_Sender < PCE_Receiver -> do all Senders see at least 1 Receiver?
# PCE_Sender == PCE_Receiver -> are all Receivers seen by both Sender and the method?
if amount_Receiver_seen_byMethod != 0:
  if len(all_sender_cells) > len(all_receiver_cells):
      ratio_Receiver_seen_byMethod = amount_Receiver_seen_byMethod / len(all_receiver_cells) * 100
  elif len(all_sender_cells) < len(all_receiver_cells):
      ratio_Receiver_seen_byMethod = amount_Receiver_seen_byMethod / len(all_sender_cells) * 100
  elif len(all_sender_cells) == len(all_receiver_cells):
      ratio_Receiver_seen_byMethod = amount_Receiver_seen_byMethod / len(all_receiver_cells) * 100
else:
  ratio_Receiver_seen_byMethod = 0

# average number of cells that each Sender cell has "seen" by the method
average_cells_perSender_seen_byMethod = (cm["weight"] != 0).sum() / len(all_sender_cells)

# ==============================================================================
# STEP 5: extract results and save
# ==============================================================================

LRdata_df = adata.uns['global_res'].sort_values("perm_pval", ascending=True)
LRdata_df = LRdata_df.loc[:, ["Ligand0", "Receptor0", "Receptor1", "perm_pval"]]
LRdata_df = LRdata_df.rename({"Ligand0": "ligand", "perm_pval": "statistics"}, axis=1)
LRdata_df = LRdata_df.sort_values("statistics")  # sort importance column in ascending order
LRdata_df["significant"] = LRdata_df["statistics"] < 0.05
LRdata_df["receptor"] = LRdata_df["Receptor0"] + "_" + LRdata_df["Receptor1"].fillna("")
LRdata_df["receptor"] = LRdata_df["receptor"].str.rstrip("_")
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]
LRdata_df["ratio_Receiver_seen_byMethod"] = ratio_Receiver_seen_byMethod
LRdata_df["average_cells_perSender_seen_byMethod"] = average_cells_perSender_seen_byMethod

LRdata_df.to_csv(significant_interactions_path, sep="\t")
