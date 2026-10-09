"""
method_mistyR.py

Runs LIANA+'s LR-MISTy (misty method, https://liana-py.readthedocs.io/en/latest/notebooks/misty.html#Ligand-Receptor-Misty)
on the semi-simulated spatial dataset and extracts the significance of the
simulated Sender -> Receiver ligand-receptor interaction.
"""

import pandas as pd
import scanpy as sc
from liana.method import lrMistyData
from liana.method.sp import LinearModel
import numpy as np
import sys
import json
import yaml
import matplotlib.pyplot as plt

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
LR_database_path = snakemake.params["LR_database"]

radius = config["radius_param"]["mistyR"][dataset][np.int64(radius_index)]

# ==============================================================================
# STEP 1: build and run the LR-MISTy model
# ==============================================================================

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

LR_database = pd.read_csv(LR_database_path, sep=" ")

misty = lrMistyData(adata,
                     resource=LR_database,
                     nz_threshold=0,   # makes the paraview also see direct neighborhoods (normally only the juxtaview would)
                     cutoff=0.01,      # doesn't influence output much
                     spatial_key='spatial',
                     kernel="gaussian",
                     bandwidth=radius,
                     use_raw=False,
                     verbose=True
                     )

misty(bypass_intra=True, model=LinearModel, verbose=True)

# ==============================================================================
# STEP 2: diagnostic plot of spatial connectivity weights, plus
# neighborhood-coverage metrics
# ==============================================================================

all_sender_cells = cm.loc[cm["Celltype_updated"] == "Sender_signalAdded", ]["Cell_ID"].tolist()
index_sender_Cells = cm.loc[cm["Celltype_updated"] == "Sender_signalAdded", ]["Cell_ID"].index.tolist()
all_receiver_cells = cm.loc[cm["Celltype_updated"] == "Receiver_signalAdded", ]["Cell_ID"].tolist()
index_receiver_cells = cm.loc[cm["Celltype_updated"] == "Receiver_signalAdded", ]["Cell_ID"].index.tolist()

# get spatial connectivities computed by MISTy for all Sender cells
cm["spatial_connectivities"] = pd.DataFrame.sparse.from_spmatrix(misty["extra"].obsp["spatial_connectivities"]).loc[index_sender_Cells, ].max(axis=0)

plt.scatter(cm["x"], cm["y"],
            c=cm['spatial_connectivities'], s=5)
plt.colorbar()
plt.scatter(cm.loc[cm["Cell_ID"].isin(all_sender_cells), "x"], cm.loc[cm["Cell_ID"].isin(all_sender_cells), "y"],
            c="#990099", s=5)
plt.scatter(cm.loc[cm["Cell_ID"].isin(all_receiver_cells), "x"], cm.loc[cm["Cell_ID"].isin(all_receiver_cells), "y"],
            c="#0000FF", s=5)
plt.xlabel("x_coord_um")
plt.ylabel("y_coord_um")
plt.title("purple - Sender cells | blue - Receiver cells | orange - cells seen by method | color bar - spatial connectivity values",
          fontsize=8)

plt.savefig(plot_neighbors_path, dpi=200)

# how many Receiver cells are seen by the method
amount_Receiver_seen_byMethod = (pd.DataFrame.sparse.from_spmatrix(misty["extra"].obsp["spatial_connectivities"]).loc[index_sender_Cells, index_receiver_cells].max(axis=0) != 0).sum()

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
average_cells_perSender_seen_byMethod = (cm["spatial_connectivities"] != 0).sum() / len(all_sender_cells)

# ==============================================================================
# STEP 3: extract results and save
# ==============================================================================

LRdata_df = misty.uns["interactions"]

# importances relate to positive/negative correlation. When ligand/receptor
# amounts are imbalanced we can get a negative correlation, which is still an
# interesting observation - so we take the absolute value for ranking purposes.
LRdata_df["importances"] = LRdata_df["importances"].abs()
LRdata_df["statistics"] = LRdata_df["importances"] 
LRdata_df = LRdata_df.sort_values("statistics", ascending=False)  # descending: higher importance first
LRdata_df["ligand_receptor"] = LRdata_df["predictor"] + "_" + LRdata_df["target"] 

# LR-MISTy models every receptor from all ligands in the resource, so its output
# contains all ligand x receptor combinations. Keep only pairs present in the resource
# so mistyR is comparable to the other methods.
db_pairs = set(LR_database["ligand"] + "_" + LR_database["receptor"])
LRdata_df = LRdata_df[LRdata_df["ligand_receptor"].isin(db_pairs)]

LRdata_df["significant"] = LRdata_df["importances"] > 2  # importance of 2 was used in the paper as a filtering criterion
LRdata_df["ratio_Receiver_seen_byMethod"] = ratio_Receiver_seen_byMethod
LRdata_df["average_cells_perSender_seen_byMethod"] = average_cells_perSender_seen_byMethod

LRdata_df.to_csv(significant_interactions_path, sep="\t")
