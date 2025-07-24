import os

configfile:"config.yaml"

rule all:
    input:
        run_processing=expand("data/processed/{dataset}/processed_counts_{dataset}.tsv" , 
            dataset=config["datasets"]),
        run_semiSimulation_inflateCounts=expand("output/{dataset}_semiSimulation_NB/inflated_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv", FC=config["semiSimulation"]["FC"], n_neighbors=config["semiSimulation"]["n_neighbors"], FC_nCells_expressingLR=config["semiSimulation"]["FC_nCells_expressingLR"],dataset=config["datasets"], indexLR_toSample=config["indexLR_toSample"]),
        run_normalization=expand("output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv" ,FC=config["semiSimulation"]["FC"], n_neighbors=config["semiSimulation"]["n_neighbors"],FC_nCells_expressingLR=config["semiSimulation"]["FC_nCells_expressingLR"], dataset=config["datasets"], indexLR_toSample=config["indexLR_toSample"]),
        run_methods=expand("output/{dataset}/{method}/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv" ,FC=config["semiSimulation"]["FC"], n_neighbors=config["semiSimulation"]["n_neighbors"],FC_nCells_expressingLR=config["semiSimulation"]["FC_nCells_expressingLR"], dataset=config["datasets"], method=config["methods"], l_param_index=config["l_param_index"], indexLR_toSample=config["indexLR_toSample"]),
        run_metrics=expand("output/{dataset}/metrics/{method}/f1_score_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.csv",FC=config["semiSimulation"]["FC"], n_neighbors=config["semiSimulation"]["n_neighbors"], FC_nCells_expressingLR=config["semiSimulation"]["FC_nCells_expressingLR"],dataset=config["datasets"], method=config["methods"], l_param_index=config["l_param_index"], indexLR_toSample=config["indexLR_toSample"])
    output:
        visualization=directory("output/results")
    shell:
        'singularity exec --no-home sing_container/liana_edgeR.sif Rscript scripts/visualization.R --config.yaml config.yaml --path_output_dir output --path_results_dir "output/results"'
            
################### Processing raw files

rule run_processing:
    threads: 1
    resources:
        mem_mb=5000
    input:
        counts="data/{dataset}/counts_{dataset}.tsv",
        metadata="data/{dataset}/metadata_{dataset}.tsv"
    output:
        processed_counts="data/processed/{dataset}/processed_counts_{dataset}.tsv",
        genemetadata="data/processed/{dataset}/genemetadata_{dataset}.RDS",
        cellmetadata="data/processed/{dataset}/cellmetadata_{dataset}.json",
        plot_allneighbors="data/processed/{dataset}/plot_allneighbors_{dataset}.pdf",
        diagnostic_plots="data/processed/{dataset}/diagnostic_plots_{dataset}.pdf"
    params:
      LR_database=config["LR_database"],
      n_neighbors=config["semiSimulation"]["n_neighbors"]
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/processing_dataset.R"
        
################### Semi-simulation
rule run_semiSimulation_inflateCounts:
    threads: 1
    resources:
        mem_mb=5000
    input:
        processed_counts="data/processed/{dataset}/processed_counts_{dataset}.tsv",
        genemetadata="data/processed/{dataset}/genemetadata_{dataset}.RDS",
        cellmetadata="data/processed/{dataset}/cellmetadata_{dataset}.json"
    output:
        inflated_counts="output/{dataset}_semiSimulation_NB/inflated_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        simulated_cellmetadata="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.json",
        simulated_interactions="output/{dataset}_semiSimulation_NB/simulated_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.RDS",
        FC_after_simulation="output/{dataset}_semiSimulation_NB/realFC_aftersimulation_{FC}_n_neigbors_{n_neighbors}_indexLR_FC_nCells_expressingLR_{FC_nCells_expressingLR}_{indexLR_toSample}.RDS",
        metadata_cpdbv5="data/cpdbv5_extrafiles/{dataset}/metadata_cpdbv5_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="data/processed/{dataset}/plot_selected_neighbors_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.pdf"
    params:
        indexLR_toSample=config["indexLR_toSample"],
    	LR_database=config["LR_database"]
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/processing_semisimulation_nbsampling.R"

################### Normalize inflated counts

rule run_normalization:
    threads: 1
    resources:
        mem_mb=5000
    input:
        inflated_counts="output/{dataset}_semiSimulation_NB/inflated_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv"
    output:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv"
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/processing_normalization.R"
        
################### Run cell-cell communication methods

rule run_method_cellphonedb:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        microenvironment="data/cpdbv5_extrafiles/microenvironment.tsv",
        metadata="data/cpdbv5_extrafiles/{dataset}/metadata_cpdbv5_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        cpdb_database="data/cpdbv5_extrafiles/cellphonedb.zip",
    output:
        significant_interactions="output/{dataset}/cellphonedb/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv"
    container:
        "sing_container/cellphonedbv5.sif"
    script:
        "scripts/method_cellphonedb.py"
        
rule run_method_lianaP_morans:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{dataset}/lianaP_morans/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{dataset}/lianaP_morans/plot_selected_neighbors_lianaP_morans_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/lianaPlus.sif"
    params:
      l_index = "{l_param_index}",
      dataset = "{dataset}",
      LR_database=config["LR_database"]
    script:
        "scripts/method_lianaP_morans.py"
        
rule run_method_lianaP_lee:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{dataset}/lianaP_lee/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{dataset}/lianaP_lee/plot_selected_neighbors_lianaP_lee_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/lianaPlus.sif"
    params:
      l_index = "{l_param_index}",
      dataset = "{dataset}",
      LR_database=config["LR_database"]
    script:
        "scripts/method_lianaP_lee.py"
        
rule run_method_spatialdm:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{dataset}/spatialdm/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{dataset}/spatialdm/plot_selected_neighbors_spatialdm_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/spatialdm.sif"
    params:
      l_index = "{l_param_index}",
      dataset = "{dataset}",
      LR_database=config["LR_database"]
    script:
        "scripts/method_spatialdm.py"
        
        
rule run_method_mistyR:
    threads: 2
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{dataset}/mistyR/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{dataset}/mistyR/plot_selected_neighbors_mistyR_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/lianaPlus.sif"
    params:
      l_index = "{l_param_index}",
      dataset = "{dataset}",
      LR_database=config["LR_database"]
    script:
        "scripts/method_mistyR.py"
        
rule run_method_NICHES:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{dataset}/NICHES/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{dataset}/NICHES/plot_selected_neighbors_NICHES_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/NICHES.sif"
    params:
      radius_index = "{l_param_index}", # NICHES uses euclidean radius for neighbor estimation
      dataset = "{dataset}",
      LR_database=config["LR_database"]
    script:
        "scripts/method_NICHES.R"
        
rule run_method_cellchat:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.json"
    output:
        significant_interactions="output/{dataset}/cellchat/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv",
        plot_neighbors="output/{dataset}/cellchat/plot_selected_neighbors_cellchat_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.png"
    container:
        "sing_container/cellchat.sif"
    params:
      radius_index = "{l_param_index}", # cellchat uses euclidean radius for neighbor estimation
      dataset = "{dataset}",
      LR_database=config["LR_database"]
    script:
        "scripts/method_cellchat.R"
################## METRICS
# f1score
rule run_metric_f1score_rankingLRgenes:
    threads: 1
    resources:
        mem_mb=1000
    input:
        significant_interactions="output/{dataset}/{method}/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv",
        simulated_interactions="output/{dataset}_semiSimulation_NB/simulated_interactions_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_indexLR_{indexLR_toSample}.RDS"
    output:
    	ranking_LRgenes="output/{dataset}/metrics/{method}/ranking_LRgenes_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.tsv",
        CT_statistics="output/{dataset}/metrics/{method}/f1_score_FC_{FC}_n_neigbors_{n_neighbors}_FC_nCells_expressingLR_{FC_nCells_expressingLR}_l_{l_param_index}_indexLR_{indexLR_toSample}.csv"
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/metric_f1score_rankingLRgenes.R"
