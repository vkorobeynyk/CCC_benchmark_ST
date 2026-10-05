import os

# ============================================================================
# Snakefile - spatial transcriptomics CCC (cell-cell communication) benchmark
#
# Pipeline overview (one pass through the DAG, per dataset x strategy):
#   1. run_processing                    - QC, Sender/Receiver cell selection,
#                                           neighbor pairing, edgeR mean/disp
#                                           estimation
#                                           (processing_dataset_spatialScattering.R
#                                           or processing_dataset_spatialColocalization.R)
#   2. run_semiSimulation_inflateCounts  - inflate ligand/receptor expression
#                                           in Sender/Receiver cells
#                                           (processing_semisimulation_nbsampling.R)
#   3. run_normalization                 - log-normalize the inflated counts
#   4. run_method_*                      - run each CCC method on the
#                                           normalized, semi-simulated data
#   5. run_metric_f1score_rankingLRgenes_all - score every method/parameter
#                                           combination against the known
#                                           simulated ground truth
#   6. visualize_results                 - generate all summary figures
#
# The wildcards {PCE_Sender}, {PCE_Receiver}, {indexLR_toSample}, and
# {radius_param_index} thread the semi-simulation and neighborhood parameters
# through every rule (see config.yaml for their meaning and value grids).
# ============================================================================

configfile: "config.yaml"

rule all:
    input:
        expand("output/{strategy}/visualize_results.done", strategy=config["strategies"])


###############################################################################
# STEP 1: Processing raw files
#
# QC + Sender/Receiver cell selection + neighbor pairing + edgeR mean/dispersion
# estimation, per dataset x strategy. The `script:` path below resolves the
# {strategy} wildcard directly, so this one rule runs either
# processing_dataset_spatialScattering.R or processing_dataset_spatialColocalization.R
# depending on which strategy Snakemake is currently building.
###############################################################################

rule run_processing:
    threads: 1
    resources:
        mem_mb=25000
    input:
        counts="data/{dataset}/counts_{dataset}.tsv",
        metadata="data/{dataset}/metadata_{dataset}.tsv"
    output:
        processed_counts="data/processed/{strategy}/{dataset}/processed_counts_{dataset}.tsv",
        genemetadata="data/processed/{strategy}/{dataset}/genemetadata_{dataset}.RDS",
        cellmetadata="data/processed/{strategy}/{dataset}/cellmetadata_{dataset}.json",
        plot_allneighbors="data/processed/{strategy}/{dataset}/plot_allneighbors_{dataset}.pdf",
        plot_radius="data/processed/{strategy}/{dataset}/plot_allradius_{dataset}.png",
        plot_SenderReceiverDistance="data/processed/{strategy}/{dataset}/plot_SenderReceiverDistance_{dataset}.pdf",
        diagnostic_plots="data/processed/{strategy}/{dataset}/diagnostic_plots_{dataset}.pdf"
    params:
        LR_database=config["LR_database"],
        dataset="{dataset}"
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/processing_dataset_{wildcards.strategy}.R"

###############################################################################
# STEP 2: Semi-simulation
#
# Inflates the expression of one ligand-receptor pair (indexLR_toSample) in a
# subset of Sender cells (ligand) and Receiver cells (receptor), sized
# according to PCE_Sender / PCE_Receiver. Also writes the CellPhoneDB-specific
# metadata format needed by run_method_cellphonedbv5.
###############################################################################

rule run_semiSimulation_inflateCounts:
    threads: 1
    resources:
        mem_mb=15000
    input:
        processed_counts="data/processed/{strategy}/{dataset}/processed_counts_{dataset}.tsv",
        genemetadata="data/processed/{strategy}/{dataset}/genemetadata_{dataset}.RDS",
        cellmetadata="data/processed/{strategy}/{dataset}/cellmetadata_{dataset}.json"
    output:
        inflated_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        simulated_cellmetadata="output/{strategy}/{dataset}_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.json",
        simulated_interactions="output/{strategy}/{dataset}_semiSimulation_NB/simulated_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.RDS",
        metadata_cpdbv5="data/cpdbv5_extrafiles/{strategy}/{dataset}/metadata_cpdbv5_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="data/processed/{strategy}/{dataset}/plot_selected_neighbors_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.png"
    params:
        LR_database=config["LR_database"]
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/processing_semisimulation_nbsampling.R"

###############################################################################
# STEP 3: Normalize inflated counts
#
# Log-normalizes the semi-simulated counts (scater::logNormCounts) so every
# CCC method downstream consumes the same normalized expression matrix.
###############################################################################

rule run_normalization:
    threads: 1
    resources:
        mem_mb=15000
    input:
        inflated_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv"
    output:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv"
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/processing_normalization.R"

###############################################################################
# STEP 4: Run cell-cell communication methods
#
# Each rule below runs one CCC method on the same normalized, semi-simulated
# input and writes out a `significant_interactions` table (ligand_receptor,
# significant, statistics, plus bookkeeping columns on how many Receiver
# cells / how much of the neighborhood the method could actually "see").
###############################################################################

rule run_method_cellphonedbv5:
    threads: 1
    resources:
        mem_mb=10000
    input:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        microenvironment="data/cpdbv5_extrafiles/microenvironment.tsv",
        cellmetadata_post_simulation="data/cpdbv5_extrafiles/{strategy}/{dataset}/metadata_cpdbv5_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        cpdb_database="data/cpdbv5_extrafiles/cellphonedb_08_25_2025_151500.zip",
    output:
        significant_interactions="output/{strategy}/{dataset}/cellphonedbv5/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv"
    container:
        "sing_container/cellphonedbv5.sif"
    script:
        "scripts/method_cellphonedbv5.py"

rule run_method_lianaP_morans:
    threads: 1
    resources:
        mem_mb=10000
    input:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{strategy}/{dataset}_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{strategy}/{dataset}/lianaP_morans/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{strategy}/{dataset}/lianaP_morans/plot_selected_neighbors_lianaP_morans_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/lianaPlus.sif"
    params:
        PCE_Receiver = "{PCE_Receiver}",
        PCE_Sender = "{PCE_Sender}",
        indexLR = "{indexLR_toSample}",
        radius_index = "{radius_param_index}",
        dataset = "{dataset}",
        LR_database=config["LR_database"]
    script:
        "scripts/method_lianaP_morans.py"

rule run_method_spatialdm:
    threads: 1
    resources:
        mem_mb=10000
    input:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{strategy}/{dataset}_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{strategy}/{dataset}/spatialdm/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{strategy}/{dataset}/spatialdm/plot_selected_neighbors_spatialdm_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/spatialdm.sif"
    params:
        radius_index = "{radius_param_index}",
        dataset = "{dataset}",
        LR_database=config["LR_database"]
    script:
        "scripts/method_spatialdm.py"


rule run_method_mistyR:
    threads: 2
    resources:
        mem_mb=10000
    input:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{strategy}/{dataset}_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{strategy}/{dataset}/mistyR/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{strategy}/{dataset}/mistyR/plot_selected_neighbors_mistyR_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/lianaPlus.sif"
    params:
        radius_index = "{radius_param_index}",
        dataset = "{dataset}",
        LR_database=config["LR_database"]
    script:
        "scripts/method_mistyR.py"

rule run_method_NICHES:
    threads: 1
    resources:
        mem_mb=35000
    input:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{strategy}/{dataset}_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{strategy}/{dataset}/NICHES/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{strategy}/{dataset}/NICHES/plot_selected_neighbors_NICHES_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/NICHES.sif"
    params:
        radius_index = "{radius_param_index}",
        dataset = "{dataset}",
        LR_database=config["LR_database"]
    script:
        "scripts/method_NICHES.R"

rule run_method_cellchat:
    threads: 1
    resources:
        mem_mb=10000
    input:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{strategy}/{dataset}_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{strategy}/{dataset}/cellchat/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{strategy}/{dataset}/cellchat/plot_selected_neighbors_cellchat_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/cellchat.sif"
    params:
        radius_index = "{radius_param_index}", # cellchat uses euclidean radius for neighbor estimation
        dataset = "{dataset}",
        LR_database=config["LR_database"]
    script:
        "scripts/method_cellchat.R"

# rule run_method_stlearn is invoked via `shell:` rather than Snakemake's usual
# `script:` directive, because of an environment/container incompatibility
# between Snakemake and stLearn (which requires python 3.10). Using shell
# directly with an argparse CLI (see method_stlearn.py) was the only way found
# to get stLearn to run inside this pipeline.
rule run_method_stlearn:
    threads: 4 
    resources:
        mem_mb=10000
    input:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{strategy}/{dataset}_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{strategy}/{dataset}/stlearn/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{strategy}/{dataset}/stlearn/plot_selected_neighbors_stlearn_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.png"
    params:
        radius_index = "{radius_param_index}",
        dataset = "{dataset}",
        LR_database=config["LR_database"]
    shell:
        "singularity exec sing_container/stlearn.sif python scripts/method_stlearn.py --normalized_counts {input.normalized_counts} --cellmetadata_post_simulation {input.cellmetadata_post_simulation} --significant_interactions {output.significant_interactions} --plot_neighbors {output.plot_neighbors} --radius_index {params.radius_index} --dataset {params.dataset} --LR_database {params.LR_database} --n_cpus {threads}"

# Baseline method: a plain Seurat Wilcoxon test with no spatial-distance
# awareness at all, to gauge how much the spatial methods actually gain from
# using spatial information.
rule run_method_seurat_wilcoxon:
    threads: 1
    resources:
        mem_mb=10000
    input:
        normalized_counts="output/{strategy}/{dataset}_semiSimulation_NB/inflated_normalized_counts_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{strategy}/{dataset}_semiSimulation_NB/simulated_cellmetadata_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{strategy}/{dataset}/seurat_wilcoxon/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv"
    container:
        "sing_container/liana_edgeR.sif"
    params:
        radius_index = "{radius_param_index}",
        dataset = "{dataset}",
        LR_database=config["LR_database"]
    script:
        "scripts/method_seurat_wilcoxon.R"

###############################################################################
# STEP 5: Metrics
#
# f1score: pulls in every (dataset, method, PCE_Sender, PCE_Receiver,
# radius_param_index, indexLR_toSample) combination's significant_interactions
# file *for one strategy*, checks whether the known simulated LR pair was
# recovered, and computes precision/recall/F1/rank for each combination.
###############################################################################

rule run_metric_f1score_rankingLRgenes_all:
    threads: 1 # workflow.cores
    resources:
        mem_mb=2000
    input:
        significant_interactions=lambda wildcards: expand(
            "output/{strategy}/{dataset}/{method}/significant_interactions_PCE_Sender_{PCE_Sender}_PCE_Receiver_{PCE_Receiver}_l_{radius_param_index}_indexLR_{indexLR_toSample}.tsv",
            strategy=wildcards.strategy,
            PCE_Sender=config["semiSimulation"]["PCE_Sender"],
            PCE_Receiver=config["semiSimulation"]["PCE_Receiver"],
            dataset=config["datasets"],
            method=config["methods"],
            radius_param_index=config["radius_param_index"],
            indexLR_toSample=config["indexLR_toSample"]
        )
    output:
        final_scores="output/{strategy}/final_scores.RDS"
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/metric_f1score_rankingLRgenes.R"

################ visualize results
rule visualize_results:
    input:
        final_scores="output/{strategy}/final_scores.RDS",
        script="scripts/visualization.R",
        config="config.yaml"
    output:
        touch("output/{strategy}/visualize_results.done")
    shell:
        """
    singularity exec sing_container/liana_edgeR.sif \
        Rscript {input.script} \
        --config.yaml {input.config} \
        --strategy {wildcards.strategy}
        """
