import pandas as pd
import scanpy as sc
import liana as li
import numpy as np
import sys
import json
import yaml
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
indexLR = snakemake.params["indexLR"]
FC_nReceiverCells = snakemake.params["FC_nReceiverCells"]
FC_nSenderCells = snakemake.params["FC_nSenderCells"]
LR_database_path = snakemake.params["LR_database"]

l = config["l_param"]["moransI"][dataset][np.int64(l_index)]
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
          fontsize = 8)
            
plt.savefig(plot_neighbors_path, dpi = 200) 

# how many CT2 cells are seen by method
amount_CT2_seen_byMethod = (cm[cm["Cell_ID"].isin(all_receiver_cells)]["spatial_connectivities"]!= 0).sum() 

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
average_cells_perCT1_seen_byMethod = (cm["spatial_connectivities"]!=0).sum() / len(all_sender_cells)

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
LRdata_df = LRdata_df[LRdata_df["morans"] > 0] # select only positive morans values as this means co-occurence
LRdata_df = LRdata_df.loc[:  , ["ligand","receptor","morans_pvals"]]
LRdata_df["significant"] = LRdata_df["morans_pvals"] < 0.05
LRdata_df = LRdata_df.rename({"morans_pvals":"statistics"},axis=1)
LRdata_df = LRdata_df.sort_values("statistics") # sort importance column on ascending order
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]
LRdata_df["ratio_CT2_seen_byMethod"] = ratio_CT2_seen_byMethod
LRdata_df["average_cells_perCT1_seen_byMethod"] = average_cells_perCT1_seen_byMethod

LRdata_df.to_csv(significant_interactions_path, sep = "\t")

############################################################
### generate plots which candidates moransI detects well ###
############################################################
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

    ax.invert_yaxis()  # common for spatial transcriptomics
    ax.set_title(f"{gene} expression (n={mask.sum()})")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(f"{gene} expression")

    plt.tight_layout()
    # DO NOT call plt.show() here
    return fig

# select 2 genes on the bottom and 2 genes on top of list
n = len(LRdata_df) -1
genes = [LRdata_df["ligand_receptor"][n].split("_")[0], LRdata_df["ligand_receptor"][n].split("_")[1], 
         LRdata_df["ligand_receptor"][0].split("_")[0], LRdata_df["ligand_receptor"][0].split("_")[1]]
if FC_nReceiverCells == 0.4 and FC_nSenderCells == 0.4 and indexLR == 1:
  with PdfPages("output/" + dataset + "/lianaP_morans/plot_spatialDistribution_LR.pdf") as pdf:
    for gene in genes:
        # Create the plot
        plot_gene_spatial(adata, gene)
        pdf.savefig() 
        plt.close()
