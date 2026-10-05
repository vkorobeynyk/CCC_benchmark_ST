"""
method_lianaP_morans.py

Runs LIANA+'s bivariate spatial Moran's I statistic on the semi-simulated
spatial dataset and extracts the significance of the simulated
Sender -> Receiver ligand-receptor interaction.
"""

import pandas as pd
import scanpy as sc
import liana as li
import numpy as np
import sys
import json
import yaml
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
indexLR = snakemake.params["indexLR"]
PCE_Receiver = snakemake.params["PCE_Receiver"]
PCE_Sender = snakemake.params["PCE_Sender"]
LR_database_path = snakemake.params["LR_database"]

radius = config["radius_param"]["moransI"][dataset][np.int64(radius_index)]

# ==============================================================================
# STEP 1: build the spatial neighbor graph
# ==============================================================================

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# build the spatial graph with the chosen kernel bandwidth
li.ut.spatial_neighbors(adata, bandwidth=radius, kernel='gaussian', set_diag=True)

# ==============================================================================
# STEP 2: diagnostic plot of spatial connectivity weights, plus
# neighborhood-coverage metrics
# ==============================================================================

all_sender_cells = cm.loc[cm["Celltype_updated"] == "Sender_signalAdded", ]["Cell_ID"].tolist()
index_sender_Cells = cm.loc[cm["Celltype_updated"] == "Sender_signalAdded", ]["Cell_ID"].index.tolist()
all_receiver_cells = cm.loc[cm["Celltype_updated"] == "Receiver_signalAdded", ]["Cell_ID"].tolist()
index_receiver_cells = cm.loc[cm["Celltype_updated"] == "Receiver_signalAdded", ]["Cell_ID"].index.tolist()

# get spatial connectivities
cm["spatial_connectivities"] = pd.DataFrame.sparse.from_spmatrix(adata.obsp["spatial_connectivities"]).loc[index_sender_Cells, ].max(axis=0)

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
amount_Receiver_seen_byMethod = (cm[cm["Cell_ID"].isin(all_receiver_cells)]["spatial_connectivities"] != 0).sum()

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
# STEP 3: run the bivariate Moran's I test and save results
# ==============================================================================

LR_database = pd.read_csv(LR_database_path, sep=" ")

liana = li.mt.bivariate(adata,
                         resource=LR_database,
                         local_name=None,    # local metric - we dont use for now as it only retrieves mean/std values
                         global_name="morans",   # global metric name
                         n_perms=200,            # number of permutations used to calculate a p-value
                         mask_negatives=False,   # whether to mask LowLow/NegativeNegative interactions
                         add_categories=True,    # whether to add local categories to the results
                         nz_prop=0,              # minimum expression proportion for ligands/receptors and their subunits
                         seed=1,
                         use_raw=False,
                         verbose=True
                         )

LRdata_df = liana
LRdata_df = LRdata_df[LRdata_df["morans"] > 0]  # keep only positive Moran's I values (co-occurrence)
LRdata_df = LRdata_df.loc[:, ["ligand", "receptor", "morans_pvals"]]
LRdata_df["significant"] = LRdata_df["morans_pvals"] < 0.05
LRdata_df = LRdata_df.rename({"morans_pvals": "statistics"}, axis=1)
LRdata_df = LRdata_df.sort_values("statistics")  # sort importance column in ascending order
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]
LRdata_df["ratio_Receiver_seen_byMethod"] = ratio_Receiver_seen_byMethod
LRdata_df["average_cells_perSender_seen_byMethod"] = average_cells_perSender_seen_byMethod

LRdata_df.to_csv(significant_interactions_path, sep="\t")

# ==============================================================================
# STEP 4: (one reference parameter combination only) plot the spatial
# distribution of the top/bottom scoring genes, to visually sanity-check
# which candidates Moran's I detects well
# ==============================================================================
'''
def plot_gene_spatial(adata, gene, min_expr=0.0, figsize=(6, 6)):
    # Check gene exists
    if gene not in adata.var_names:
        raise ValueError(f"Gene '{gene}' not found in adata.var_names")

    # Extract spatial coordinates
    coords = adata.obsm["spatial"]
    x = coords[:, 0]
    y = coords[:, 1]

    # Extract gene expression
    expr = adata[:, gene].X.A.flatten() if hasattr(adata[:, gene].X, "A") else np.array(adata[:, gene].X).flatten()

    # Filter cells with expression > min_expr
    mask = expr > min_expr

    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # Plot all cells in gray
    ax.scatter(x, y, s=5, color="lightgray", alpha=0.4)

    # Plot expressing cells colored by expression intensity
    sc = ax.scatter(
        x[mask],
        y[mask],
        c=expr[mask],
        s=10,
        cmap="viridis",
        edgecolors="none"
    )

    ax.invert_yaxis()  # common convention for spatial transcriptomics
    ax.set_title(f"{gene} expression (n={mask.sum()})")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(f"{gene} expression")

    plt.tight_layout()
    # DO NOT call plt.show() here
    return fig

# select the bottom 2 and top 2 genes from the ranked list
n = len(LRdata_df) - 1
genes = [LRdata_df["ligand_receptor"][n].split("_")[0], LRdata_df["ligand_receptor"][n].split("_")[1],
         LRdata_df["ligand_receptor"][0].split("_")[0], LRdata_df["ligand_receptor"][0].split("_")[1]]
if PCE_Receiver == 2 and PCE_Sender == 2 and indexLR == 1:
    with PdfPages("output/" + dataset + "/lianaP_morans/plot_spatialDistribution_LR.pdf") as pdf:
        for gene in genes:
            # Create the plot
            plot_gene_spatial(adata, gene)
            pdf.savefig()
            plt.close()
'''
