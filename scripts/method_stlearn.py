"""
Runs stLearn's cell-cell communication analysis on the semi-simulated
spatial dataset and extracts the significance of the simulated
Sender -> Receiver ligand-receptor interaction.

This script is invoked directly via `shell:` in the Snakefile (rather than
Snakemake's `script:` directive) due to environment/container incompatibilities
with running stLearn (python 3.10) through the usual Snakemake R/Python hook -
hence the argparse-based CLI interface below.

CLI arguments: --normalized_counts, --cellmetadata_post_simulation,
               --significant_interactions, --plot_neighbors,
               --radius_index, --dataset, --LR_database
"""

import matplotlib.pyplot as plt
import stlearn as st
import scanpy as sc
import pandas as pd
import json
import numpy as np
import yaml
import itertools
import argparse

# ==============================================================================
# STEP 0: CLI args / Snakemake I/O and params
# ==============================================================================

parser = argparse.ArgumentParser(
    description="Run stLearn cell-cell communication analysis"
)
parser.add_argument("--normalized_counts", required=True)
parser.add_argument("--cellmetadata_post_simulation", required=True)
parser.add_argument("--significant_interactions", required=True)
parser.add_argument("--plot_neighbors", required=True)
parser.add_argument("--radius_index", required=True)
parser.add_argument("--dataset", required=True)
parser.add_argument("--LR_database", required=True)
parser.add_argument("--n_cpus", type=int, default=1,
                     help="CPUs to hand to stLearn's permutation test (st.tl.cci.run). ")
args = parser.parse_args()

# input files
normalized_counts_path = args.normalized_counts
cellmetadata_path = args.cellmetadata_post_simulation

adata = sc.AnnData(pd.read_csv(normalized_counts_path, sep="\t").T)
with open(cellmetadata_path, "r") as f:
    cellmetadata = json.load(f)

cm = pd.DataFrame(cellmetadata["metadata"])

# output files
significant_interactions_path = args.significant_interactions
plot_neighbors_path = args.plot_neighbors

# params
with open("config.yaml", "r") as stream:
    config = yaml.safe_load(stream)

radius_index = args.radius_index
dataset = args.dataset
LR_database_path = args.LR_database

radius = config["radius_param"]["stlearn"][dataset][np.int64(radius_index)]

adata.obs[["imagerow", "imagecol"]] = np.array([cm["x"], cm["y"]]).T
adata.obs["Celltype"] = list(cm["Celltype"])

# ==============================================================================
# STEP 1: expand the LR database to single ligand-receptor gene pairs
#
# stlearn doesn't account for subunits and requires plain L-R nomenclature.
# Since other methods use subunit information, we split L-R1-R2 entries into
# L-R1 and L-R2 here and aggregate the results back together at the end.
# ==============================================================================

lrs = pd.read_csv(LR_database_path, sep=" ")
s = lrs["ligand_receptor"]

def expand_entry(entry):
    parts = entry.split("_")
    if len(parts) == 2:
        return [(entry, entry)]
    else:
        # pair the ligand with each receptor subunit
        return [(f"{parts[0]}_{p}", entry) for p in parts[1:]]

LR_expanded = s.apply(expand_entry).explode().reset_index(drop=True)
LR_expanded = pd.DataFrame(LR_expanded.tolist(), columns=["expanded", "ligand_receptor"])
# keep only unique "expanded" values
LR_expanded = LR_expanded.drop_duplicates(subset="expanded").reset_index(drop=True)

# ==============================================================================
# STEP 2: run stLearn's gene-level and celltype-level CCI tests
# ==============================================================================

# Gene-level permutation test: identifies which LR pairs are spatially
# co-expressed more than expected by chance.
st.tl.cci.run(adata, LR_expanded["expanded"].to_numpy(),
              min_spots=0,
              distance=radius,
              n_pairs=100,
              n_cpus=args.n_cpus
              )

st.tl.cci.adj_pvals(adata, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')

# Celltype enrichment and permutation: summarizes LR pairs between annotated
# celltypes and tests whether those connections occur more often than expected.
st.tl.cci.run_cci(adata, 'Celltype',
                   min_spots=0,
                   spot_mixtures=False,

                   cell_prop_cutoff=0,
                   sig_spots=True,
                   n_perms=10,
                   n_cpus=args.n_cpus
                   )

# ==============================================================================
# STEP 3: extract Sender -> Receiver p-values and re-aggregate subunits
# ==============================================================================

# Get p-values for every interaction, for Sender -> Receiver
p_vals = []
for LR in adata.uns["per_lr_cci_pvals_Celltype"].keys():
    p_vals.append(adata.uns["per_lr_cci_pvals_Celltype"][LR].loc["Sender", "Receiver"])

cci_pvals = pd.DataFrame(dict(ligand_receptor=adata.uns["per_lr_cci_pvals_Celltype"].keys(),
                               statistics=p_vals,
                               significant=[i < 0.05 for i in p_vals]))

# Aggregate results from L-R1 | L-R2 back into L-R1-R2
rows = []
for entry in LR_expanded["ligand_receptor"]:
    # generate all possible pairwise combinations for L-R1-R2
    x = entry.split("_")
    pairs_joined = [f"{a}_{b}" for a, b in itertools.combinations(x, 2)]
    # locate the combinations in the results dataframe
    tmp_df = cci_pvals[cci_pvals.isin(pairs_joined).any(axis=1)]

    # average the statistic column over the subunit-level results
    if tmp_df.shape[0] != 0:
        rows.append({
            "ligand_receptor": entry,
            "statistics": tmp_df["statistics"].mean(),
            "significant": all(tmp_df["significant"])})

LRdata_df = pd.DataFrame(rows).drop_duplicates("ligand_receptor").reset_index(drop=True)
LRdata_df = LRdata_df.sort_values("statistics")  # sort importance column in ascending order

# ==============================================================================
# STEP 4: diagnostic plot of the method's neighborhood, plus
# neighborhood-coverage metrics
# ==============================================================================

Sender = cm.loc[cm["Celltype"] == "Sender", ]["Cell_ID"].tolist()
Receiver = cm.loc[cm["Celltype"] == "Receiver", ]["Cell_ID"].tolist()
all_sender_cells = cm.loc[cm["Celltype_updated"] == "Sender_signalAdded", ]["Cell_ID"].tolist()
all_receiver_cells = cm.loc[cm["Celltype_updated"] == "Receiver_signalAdded", ]["Cell_ID"].tolist()

# get all cells that the method "sees"
cells_seen_by_method = [item for sublist in adata.obsm["spot_neigh_bcs"]["neighbour_bcs"][all_sender_cells] for item in sublist.split(",") if item != ""]
cells_seen_by_method = list(set(cells_seen_by_method) - set(Sender))  # exclude the Sender cells themselves

plt.scatter(cm.loc[:, "x"], cm.loc[:, "y"],
            c="grey", s=2)
plt.scatter(cm.loc[cm["Cell_ID"].isin(Sender), "x"], cm.loc[cm["Cell_ID"].isin(Sender), "y"],
            c="#FFCCFF", s=8)
plt.scatter(cm.loc[cm["Cell_ID"].isin(cells_seen_by_method), "x"], cm.loc[cm["Cell_ID"].isin(cells_seen_by_method), "y"],
            c="orange", s=8)
plt.scatter(cm.loc[cm["Cell_ID"].isin(all_sender_cells), "x"], cm.loc[cm["Cell_ID"].isin(all_sender_cells), "y"],
            c="#990099", s=12)
plt.scatter(cm.loc[cm["Cell_ID"].isin(all_receiver_cells), "x"], cm.loc[cm["Cell_ID"].isin(all_receiver_cells), "y"],
            c="#0000FF", s=12)
plt.xlabel("x_coord_um")
plt.ylabel("y_coord_um")
plt.title("purple - Sender cells | blue - Receiver cells | orange - cells seen by method",
          fontsize=8)

plt.savefig(plot_neighbors_path, dpi=200)

# how many Receiver cells are seen by the method
amount_Receiver_seen_byMethod = len(set(cells_seen_by_method).intersection(set(all_receiver_cells)))

# PCE_Sender > PCE_Receiver -> how many Receivers are seen by the method?
# PCE_Sender < PCE_Receiver -> do all Senders see at least 1 Receiver?
# PCE_Sender == PCE_Receiver -> are all Receivers seen by both Sender and the method?
if amount_Receiver_seen_byMethod != 0:
    if len(all_sender_cells) > len(all_receiver_cells):
        LRdata_df["ratio_Receiver_seen_byMethod"] = amount_Receiver_seen_byMethod / len(all_receiver_cells) * 100
    elif len(all_sender_cells) < len(all_receiver_cells):
        LRdata_df["ratio_Receiver_seen_byMethod"] = amount_Receiver_seen_byMethod / len(all_sender_cells) * 100
    elif len(all_sender_cells) == len(all_receiver_cells):
        LRdata_df["ratio_Receiver_seen_byMethod"] = amount_Receiver_seen_byMethod / len(all_receiver_cells) * 100
else:
    LRdata_df["ratio_Receiver_seen_byMethod"] = 0

# average number of cells that each Sender cell has "seen" by the method
LRdata_df["average_cells_perSender_seen_byMethod"] = len(cells_seen_by_method) / len(all_sender_cells)

# ==============================================================================
# STEP 5: save results
# ==============================================================================

LRdata_df.to_csv(significant_interactions_path, sep="\t")
