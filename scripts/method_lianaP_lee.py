import pandas as pd
import scanpy as sc
import liana as li
import numpy as np
import sys
import json
from plotnine import ggplot, geom_point, aes , ggtitle

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

##############
### Params ###
##############
l = cellmetadata["spatialWeight_lianaPlus"]

#adata = sc.AnnData(pd.read_csv("output/STARmap_plus_HPC_semiSimulation_NB/inflated_normalized_counts_FC_1_n_neigbors_4.tsv", sep="\t").T)
#cellmetadata_path = "output/STARmap_plus_HPC_semiSimulation_NB/simulated_cellmetadata_1_n_neigbors_4.json"

adata.obsm["spatial"] = np.array([cm["x"], cm["y"]]).T

# build the spatial graph with a selected bandwidth
li.ut.spatial_neighbors(adata, bandwidth=l[0], kernel='gaussian', set_diag=True)

# Bivariate Ligand-Receptor Relationships
lrdata = li.mt.bivariate(adata,
                      resource_name='cellchatdb', # NOTE: uses HUMAN gene symbols!
                      local_name='cosine', # Name of the function - currenty the other local metrics dont work/ dont change result at all
                      global_name="lee", # Name global function
                      n_perms=100, # Number of permutations to calculate a p-value
                      mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                      add_categories=True, # Whether to add local categories to the results
                      nz_prop=0, # Minimum expr. proportion for ligands/receptors and their subunits
                      seed=1,
                      use_raw=False,
                      verbose=True
)

#LRadata = out[1] # subset of adata to only LR only 
LRdata_df = lrdata.var.sort_values("lee_pvals", ascending=True) # extract df with interactions and statistics

# Plot interactions
#sc.pl.spatial(LRadata ,color=['LGI3^ADAM23'], size=1.4, vmax=1, spot_size=50 ,cmap='magma')
#sc.pl.spatial(adata, color=['LGI3', 'ADAM23'], size=1.4, ncols = 2, spot_size=50)

# save data
LRdata_df = LRdata_df.loc[LRdata_df["lee_pvals"] < 0.05 , ["ligand","receptor","lee_pvals"]]
LRdata_df = LRdata_df.rename({"lee_pvals":"pval"},axis=1)
LRdata_df["ligand_receptor"] = LRdata_df["ligand"] + "_" + LRdata_df["receptor"]
LRdata_df.to_csv(significant_interactions_path, sep = "\t")
