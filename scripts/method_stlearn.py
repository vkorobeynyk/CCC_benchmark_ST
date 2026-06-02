import matplotlib.pyplot as plt
import stlearn as st
import scanpy as sc
import pandas as pd
import json
import numpy as np
import yaml
import itertools
import argparse

#⃣ Create a parser
parser = argparse.ArgumentParser(
    description="Run stLearn cell-cell communication analysis"
)

# Define expected arguments
parser.add_argument("--normalized_counts", required=True)
parser.add_argument("--cellmetadata_post_simulation", required=True)
parser.add_argument("--significant_interactions", required=True)
parser.add_argument("--plot_neighbors", required=True)
parser.add_argument("--l_index", required=True)
parser.add_argument("--dataset", required=True)
parser.add_argument("--LR_database",required=True)

# Parse the arguments
args = parser.parse_args()

#############
### INPUT ###
#############
normalized_counts_path = args.normalized_counts
print("here")
cellmetadata_path = args.cellmetadata_post_simulation

adata = sc.AnnData(pd.read_csv(normalized_counts_path, sep="\t").T)
with open(cellmetadata_path,"r") as f:
  cellmetadata = json.load(f)
  
cm = pd.DataFrame(cellmetadata["metadata"])

##############
### OUTPUT ###
##############
significant_interactions_path = args.significant_interactions
plot_neighbors_path = args.plot_neighbors

##############
### Params ###
##############
with open("config.yaml","r") as stream:
  config = yaml.safe_load(stream)

l_index = args.l_index
dataset = args.dataset
LR_database_path = args.LR_database

l = config["l_param"]["stlearn"][dataset][np.int64(l_index)]
#adata = sc.AnnData(pd.read_csv("output/MERFISH_mColon_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_2_indexLR_1.tsv", sep="\t").T)
#cellmetadata_path = "output/Slideseq2_HPC_semiSimulation_NB/simulated_cellmetadata_FC_1_FC_nSenderCells_2_FC_nReceiverCells_1_indexLR_1.json"

adata.obs[["imagerow","imagecol"]] = np.array([cm["x"], cm["y"]]).T
adata.obs["Celltype"] = list(cm["Celltype"])


# stlearn doesnt not account for subunits and requires L-R nomenclature. As other methods use subunit information I am here splitting L-R1-R2 into L-R1 and L-R2 and then
# aggregate the results int the end

# Expand LR databse from L-R1-R2 to L-R1 and L-R2
lrs = pd.read_csv(LR_database_path, sep=" ")
s = lrs["ligand_receptor"]

def expand_entry(entry):
    parts = entry.split("_")
    if len(parts) == 2:
        # Return unchanged but keep also original
        return [(entry, entry)]
    else:
        # Pair the first with each of the rest
        return [(f"{parts[0]}_{p}", entry) for p in parts[1:]]
        
LR_expanded = s.apply(expand_entry).explode().reset_index(drop=True)
LR_expanded = pd.DataFrame(LR_expanded.tolist(), columns=["expanded", "ligand_receptor"])
# Keep only unique "expanded" values
LR_expanded = LR_expanded.drop_duplicates(subset="expanded").reset_index(drop=True)

# gene level permutation
# Identifies which LR are spatially co-expressed more than expected 
st.tl.cci.run(adata, LR_expanded["expanded"].to_numpy(),
              min_spots=0,  # Filter out any LR pairs with no scores for less than min_spots
              distance=l,  # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
              n_pairs=250,  # Number of random pairs to generate; low as example, recommend ~10,000
              n_cpus=None   # Number of CPUs for parallel. If None, detects & use all available.
              )

st.tl.cci.adj_pvals(adata, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')

# celltype celltype enrichment and permutation
# Summarizes those LR pairs between annotated cell types and check if those connections occur more often than expected
st.tl.cci.run_cci(adata, 'Celltype',  # Spot cell information either in data.obs or data.uns
                  min_spots=0,          # Minimum number of spots for LR to be tested.
                  spot_mixtures=False,   # If True will use the label transfer scores,
                                       
                  cell_prop_cutoff=0, 
                  sig_spots=True,      
                  n_perms=10            
                 )
                 

# Get pvalues for every interaction for CT1 - CT2
p_vals = []
for LR in adata.uns["per_lr_cci_pvals_Celltype"].keys():
    p_vals.append(adata.uns["per_lr_cci_pvals_Celltype"][LR].loc["CT1","CT2"])

cci_pvals = pd.DataFrame(dict(ligand_receptor = adata.uns["per_lr_cci_pvals_Celltype"].keys(), 
                  statistics = p_vals,
                  significant = [i < 0.05 for i in p_vals]))

# Aggregate results from L-R1 | L-R2 to L-R1-R2
rows = []
for entry in LR_expanded["ligand_receptor"]:
    # create all posible combinations for L-R1-R2 
    x = entry.split("_")
    pairs_joined = [f"{a}_{b}" for a, b in itertools.combinations(x, 2)]
    # locate where are the combinations in the dataframe
    tmp_df = cci_pvals[cci_pvals.isin(pairs_joined).any(axis=1)]

    # For L-R that exist in the stlearn output average over the statistic column
    if tmp_df.shape[0] != 0:
        # average results 
        rows.append({
            "ligand_receptor":entry,
            "statistics":tmp_df["statistics"].mean(),
            "significant":all(tmp_df["significant"])})

LRdata_df = pd.DataFrame(rows).drop_duplicates("ligand_receptor").reset_index(drop=True)
LRdata_df = LRdata_df.sort_values("statistics") # sort importance column on ascending order

#########################################
### Plot weights according to l param ###

CT1 = cm.loc[cm["Celltype"] == "CT1",]["Cell_ID"].tolist()
CT2 = cm.loc[cm["Celltype"] == "CT2",]["Cell_ID"].tolist()
all_sender_cells = cm.loc[cm["Celltype_updated"] == "CT1_signalAdded",]["Cell_ID"].tolist()
all_receiver_cells = cm.loc[cm["Celltype_updated"] == "CT2_signalAdded",]["Cell_ID"].tolist()
# get all cells that method "sees"
cells_seen_by_method = [item for sublist in adata.obsm["spot_neigh_bcs"]["neighbour_bcs"][all_sender_cells] for item in sublist.split(",") if item != ""]
cells_seen_by_method = list(set(cells_seen_by_method) - set(CT1)) # subtract CT1 cells

plt.scatter(cm.loc[:,"x"], cm.loc[:,"y"], 
            c= "grey", s = 2)
plt.scatter(cm.loc[cm["Cell_ID"].isin(CT1),"x"], cm.loc[cm["Cell_ID"].isin(CT1),"y"], 
            c= "#FFCCFF", s = 8)
plt.scatter(cm.loc[cm["Cell_ID"].isin(cells_seen_by_method),"x"], cm.loc[cm["Cell_ID"].isin(cells_seen_by_method),"y"], 
            c= "orange", s = 8)
plt.scatter(cm.loc[cm["Cell_ID"].isin(all_sender_cells),"x"], cm.loc[cm["Cell_ID"].isin(all_sender_cells),"y"], 
            c= "#990099", s = 12)
plt.scatter(cm.loc[cm["Cell_ID"].isin(all_receiver_cells),"x"], cm.loc[cm["Cell_ID"].isin(all_receiver_cells),"y"], 
            c= "#0000FF", s = 12)
plt.xlabel("x_coord_um")
plt.ylabel("y_coord_um")
plt.title("purple - sender cells | blue - receiver cells | orange - cells seen by method",
          fontsize = 8)
            
plt.savefig(plot_neighbors_path, dpi = 200) 

# how many CT2 cells are seen by method
amount_CT2_seen_byMethod = len(set(cells_seen_by_method).intersection(set(all_receiver_cells)))

# FCsender > FCreceiver -> how many receivers are seen by method?
# FCsender < FCreceiver -> do all senders see 1 receiver?
# FCsender == FCreceiver -> are all receiver seen by CT1 and method
if amount_CT2_seen_byMethod != 0:
    if len(all_sender_cells) > len(all_receiver_cells):
        LRdata_df["ratio_CT2_seen_byMethod"] = amount_CT2_seen_byMethod / len(all_receiver_cells) * 100
    elif len(all_sender_cells) < len(all_receiver_cells):
        LRdata_df["ratio_CT2_seen_byMethod"] = len(all_sender_cells) / amount_CT2_seen_byMethod * 100
        # if there are more than 1 receiver per sender
        #if amount_CT2_seen_byMethod > len(all_sender_cells):
        #    ratio_CT2_seen_byMethod = 100
    elif len(all_sender_cells) == len(all_receiver_cells):
        LRdata_df["ratio_CT2_seen_byMethod"] = amount_CT2_seen_byMethod / len(all_receiver_cells) * 100
else:
    LRdata_df["ratio_CT2_seen_byMethod"] = 0


# average cells that each CT1 has that are seen by method
LRdata_df["average_cells_perCT1_seen_byMethod"] = len(cells_seen_by_method) / len(all_sender_cells)
        
LRdata_df.to_csv(significant_interactions_path, sep = "\t")
