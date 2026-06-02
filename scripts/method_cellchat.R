'
1- To infer the cell state-specific communications, CellChat identifies over-expressed ligands or receptors in one cell group and then 
identifies over-expressed ligand-receptor interactions if either ligand or receptor are over-expressed.
Cellchat also provides a function to project gene expression data onto protein-protein interaction (PPI) network. 
Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence 
experimentally validated protein-protein network.
2- CellChat infers the biologically significant cell-cell communication by assigning
each interaction with a probability value and peforming a permutation test
'

suppressMessages({
  # Load package
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(CellChat)
  library(jsonlite)
  library(magrittr)
  library(purrr)
  source("scripts/helper_functions.R")
})

# An useful error if the argument is missing
if (is.null(snakemake@input[["normalized_counts"]]) | is.null(snakemake@input[["cellmetadata_post_simulation"]]) 
    | is.null(snakemake@output[["plot_neighbors"]]) | is.null(snakemake@output[["significant_interactions"]]) ){
  stop("Argument_name needs to be specified, but is missing.n", call.=FALSE)
}

# Read the argument
normalized_counts_path = snakemake@input[["normalized_counts"]]
cellmetadata_path = snakemake@input[["cellmetadata_post_simulation"]]
significant_interactions_path = snakemake@output[["significant_interactions"]]
plot_neighbors_path = snakemake@output[["plot_neighbors"]]

##############
### Params ###
##############
config = yaml::read_yaml("config.yaml")

l_index = as.integer(snakemake@params["l_index"]) +1
dataset = snakemake@params["dataset"] %>% as.character
radius = config[["l_param"]][["CellChat"]][[dataset]][l_index] %>% unlist

LR_database_path = snakemake@params[["LR_database"]]
LR_database = read.table(LR_database_path, row.names = 1)

#############
# load data #
#############

inflated_counts = read.csv(normalized_counts_path,sep="\t") %>% as.matrix
cellmetadata = read_json(path = cellmetadata_path)

'
inflated_counts = read.csv("output/Visium_HD_HPC_semiSimulation_NB//inflated_normalized_counts_FC_1_FC_nSenderCells_0.5_FC_nReceiverCells_6_indexLR_2.tsv",sep="\t") %>% as.matrix
cellmetadata = read_json(path = "output/Visium_HD_HPC_semiSimulation_NB/simulated_cellmetadata_FC_1_FC_nSenderCells_0.5_FC_nReceiverCells_6_indexLR_2.json")
LR_database = read.table("data/LR_database.tsv", row.names = 1)

'

# transform json list to individual dataframe
cellmetadata = convert_json_to_df(cellmetadata)
rownames(cellmetadata$metadata) = cellmetadata$metadata$Cell_ID

########################################################
# Compute all sender and receiver cells in the dataset #
########################################################

# CODE ADAPTED FROM NICHES GITHUB TO CALCULATE THE RADIUS FOR PLOTTING REASONS
# Compute the euclidean distance matrix
coord = cellmetadata$metadata %>% 
  select(x,y) %>% 
  mutate_at(vars(x,y) , as.numeric)

distance_mat = apply(coord, 1, function(pt)
  (sqrt(abs(pt["x"] - coord$x)^2 + abs(pt["y"] - coord$y)^2))
)

# generate a list where each index name is sender cell and it contains all cells within the radius seen by the method
CT1 = cellmetadata$metadata %>% filter(Celltype == "CT1") %>% pull(Cell_ID)
CT2 = cellmetadata$metadata %>% filter(Celltype == "CT2") %>% pull(Cell_ID)
CT1_signalAdded = cellmetadata$metadata %>% filter(Celltype_updated == "CT1_signalAdded") %>% pull(Cell_ID)
CT2_signalAdded = cellmetadata$metadata %>% filter(Celltype_updated == "CT2_signalAdded") %>% pull(Cell_ID)
vec = map(CT1_signalAdded, function(cell_OI) {
  within_radius = distance_mat[,cell_OI == colnames(distance_mat)] < radius
  return(which(within_radius))
}) %>% as.list

all_cells_seen_byMethod = colnames(distance_mat)[unlist(vec) %>% unique]
all_cells_seen_byMethod = all_cells_seen_byMethod[!all_cells_seen_byMethod %in% CT1_signalAdded] # remove CT1 cells

# plot
plt = ggplot(coord, aes(x = x ,y = y)) + 
  geom_point(size = 0.1) +
  geom_point(data=coord[all_cells_seen_byMethod,] , aes(x=x, y=y), colour="orange", size=2) +
  geom_point(data=coord[CT1,] , aes(x=x, y=y), colour="#FFCCFF", size=2) +
  geom_point(data=coord[CT1_signalAdded,] , aes(x=x, y=y), colour="#990099", size=3) +
  geom_point(data=coord[CT2_signalAdded,] , aes(x=x, y=y), colour="#0000FF", size=3) +
  ggtitle("Cellchat euclidean filtering pink -> CT1 |orange -> cells seen by method | purple -> CT1_signalAdded | blue -> CT2_signalAdded")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  theme_bw()

ggsave(plot_neighbors_path, plt, device = "png", width = 30, height = 25, units = "cm")

#########################################
# How many CT2 cells are seen by method #
#########################################
amount_CT1_signalAdded_cells = cellmetadata$metadata %>% filter(Celltype_updated == "CT1_signalAdded") %>% nrow
amount_CT2_signalAdded_cells = cellmetadata$metadata %>% filter(Celltype_updated == "CT2_signalAdded") %>% nrow
amount_CT2_seen_byMethod = which(all_cells_seen_byMethod %in% (cellmetadata$metadata %>% filter(Celltype_updated == "CT2_signalAdded") %>% pull(Cell_ID))) %>% length

# FCsender > FCreceiver -> how many receivers are seen by method?
# FCsender < FCreceiver -> do all senders see 1 receiver?
# FCsender == FCreceiver -> are all receiver seen by CT1 and method
if(amount_CT2_seen_byMethod != 0) {
  if (amount_CT1_signalAdded_cells > amount_CT2_signalAdded_cells) {
    ratio_CT2_seen_byMethod = amount_CT2_seen_byMethod / amount_CT2_signalAdded_cells * 100
  } else if (amount_CT1_signalAdded_cells < amount_CT2_signalAdded_cells) {
    ratio_CT2_seen_byMethod = amount_CT1_signalAdded_cells / amount_CT2_seen_byMethod * 100
  } else if (amount_CT1_signalAdded_cells == amount_CT2_signalAdded_cells) {
    ratio_CT2_seen_byMethod = amount_CT2_seen_byMethod / amount_CT2_signalAdded_cells * 100
  }
} else {ratio_CT2_seen_byMethod = 0}

# average cells that each CT1 has that are seen by method
average_cells_perCT1_seen_byMethod = length(all_cells_seen_byMethod) / amount_CT1_signalAdded_cells

##############
# Run method #
##############
# Cellchat, if needed, transforms data from pixels to micrometers. For that, it uses 2 parameters: ratio and tol
# Ratio is a conversion rate from pixels to micrometers. As all our data is in micrometers, ratio = 1.
# Tol is a factor that is important if one compares center-to-center distance against "interaction range" parameter because tol is usually half of cell/spot size
# as the goal is just for methods to use same distance in space, tol is not important. It is set to 1
# all datasets were rescaled to micrometers
spatial.factors = data.frame(ratio = 1, tol = 1)

# Create a CellChat object
cellchat = createCellChat(object = inflated_counts, meta = data.frame(Celltype = cellmetadata$metadata$Celltype,
                                                                   samples = "sample1" %>% as.factor, 
                                                                   row.names = cellmetadata$metadata$Cell_ID), 
                           group.by = "Celltype",datatype = "spatial", coordinates = as.matrix(coord), spatial.factors = spatial.factors)

CellChatDB = CellChatDB.human

################################################################
### Modify database of CellChat according to LR_database.tsv ###
################################################################
# Here I am filtering CellchatDB according to LR_database.tsv file
'
# To understand how CellchatDB works, I looked at the example (from LR_database.tsv) of:
# ligand (L) -> WNT4
# receptor (R) -> LRP6_FZD8
# Cellchatdb has the interaction as WNT4_FZD8_LRP6 but LR_database has WNT4_LRP6_FZD8
# Here I switch the receptor order
# I cant switch the order in cellchatDB directly because otherwise the cellchat wont recognize the interaction as valid

fix_gene_order = function(df) {
  sapply(df, function(x) 
  {
    parts = strsplit(x, "_")[[1]]
    if(length(parts) == 2) {str_c(parts[2] , "_" , parts[1])
    } else {x}
  })
}

LR_database$receptor = fix_gene_order(LR_database$receptor)
LR_database$ligand = fix_gene_order(LR_database$ligand)
LR_database$interaction_name  = str_c(LR_database$ligand, "_", LR_database$receptor)


just the same as in LR_database.tsv file
'
# cellchatDB is a list with 4 entries:
# Interactions -> I subset to same ones as LR_database.tsv
# geneInfo has information of single genes -> not necessary to filter
# complex -> no need to filter because cellchat just fetches what it needs
# cofactors I am not simulating -> I simply remove them from the database
CellChatDB$interaction = CellChatDB$interaction[which(CellChatDB$interaction$interaction_name %in% LR_database$ligand_receptor),]

stopifnot(nrow(CellChatDB$interaction) == nrow(LR_database))

# replace the cofactor dataframe by empty strings 
x = CellChatDB$cofactor %>% apply(.,2, function(x) {return(rep("", length(x)))})
rownames(x) = rownames(CellChatDB$cofactor)
CellChatDB$cofactor = x %>% as.data.frame() # cellchat requires dataframe

cellchat@DB = CellChatDB # dim() same as in LRdatabase

cellchat = subsetData(cellchat) # This step is necessary even if using the whole database

# wilcoxon test to remove features
# only uses pvalue threshold
cellchat = identifyOverExpressedGenes(cellchat,min.cells = 0,thresh.fc = 0,thresh.p = 0.05) 
cellchat = identifyOverExpressedInteractions(cellchat) # 

# if wilcoxon test results didnt find any significant LR to test for spatial
if(nrow(cellchat@LR$LRsig) == 0)
{
  write.table(data.frame(ligand_receptor = NA , significant = FALSE, statistics = 0, 
                         ratio_CT2_seen_byMethod = ratio_CT2_seen_byMethod, 
                         average_cells_perCT1_seen_byMethod = average_cells_perCT1_seen_byMethod) ,significant_interactions_path,
              sep = "\t")
} else{
  '
  When inferring contact-dependent or juxtacrine signaling, users should provide a value of contact.range and set contact.dependent = TRUE. 
  Briefly, users can set contact.range = 10, which is a typical human cell size. 
  However, for low-resolution spatial data such as 10X visium, it should be the cell center-to-center distance (i.e., contact.range = 100 for 10X visium data). 
  Please check the vignette of FAQ on applying CellChat to spatially resolved transcriptomics data for detailed explanations. 
  In this example, we did not use the L-R pairs from Cell-Cell Contact signaling, therefore we can set contact.dependent = FALSE and contact.range = NULL. 
  But as an illustration, we use the following settings that lead to the same results.
  
  Of note, ‘trimean’ approximates 25% truncated mean, 
  implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%
  '
    # When comparing communication across different CellChat objects, the same scale factor should be used
    # I tested different scale.distance (0.1,0.5,1) and the results didnt change. What changes is the probability but pvalue always stays the same
    
    '
  Re: nboot  (https://github.com/sqjin/CellChat/issues/244)
  I think the results will not change too much. If nboot = 100, then thresh = 0.05 means there are five permuations having larger 
  communication probabilities. If nboot = 20, then thresh = 0.05 means there are one permutation having larger communication pprobabilities.
  '
    
  # As we are testing amount of cells that should express each gene, I set trim = 0.001 -> 0.1% of cells have to express the gene 
  cellchat = computeCommunProb(cellchat, type = "truncatedMean", trim = 0.001,
                               distance.use = TRUE, interaction.range = radius, scale.distance = 1, 
                               contact.dependent = FALSE,contact.range = NULL, nboot = 100)
  
  df.net = subsetCommunication(cellchat, thresh = 1) # threshold of the p-value for determining significant interaction
  
  df.net %<>% filter(source == "CT1" & target == "CT2") %>% mutate(ligand_receptor = gsub("—","_",interaction_name),
                                                                   significant = pval < 0.05, # only returns significant 
                                                                   statistics = prob) %>% dplyr::arrange(desc(prob))
  
  
  # save data
  if(nrow(df.net) != 0)
  {
    write.table(data.frame(ligand_receptor = df.net$interaction_name , significant = df.net$significant, statistics = df.net$prob, 
                           ratio_CT2_seen_byMethod = ratio_CT2_seen_byMethod,
                           average_cells_perCT1_seen_byMethod = average_cells_perCT1_seen_byMethod) ,significant_interactions_path,
                sep = "\t")
  } else {
    write.table(data.frame(ligand_receptor = NA , significant = FALSE, statistics = 0, 
                           ratio_CT2_seen_byMethod = ratio_CT2_seen_byMethod, 
                           average_cells_perCT1_seen_byMethod = average_cells_perCT1_seen_byMethod) ,significant_interactions_path,
                sep = "\t")
  }
}

