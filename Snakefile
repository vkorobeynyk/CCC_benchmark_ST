import os

configfile:"config.yaml"

rule all:
    input:
        run_processing=expand("data/processed/{dataset}/processed_counts_{dataset}.tsv" , 
            dataset=config["datasets"]),
        run_semiSimulation_inflateCounts=expand("output/{dataset}_semiSimulation_NB/inflated_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv", FC=config["semiSimulation"]["FC"], n_neighbors=config["semiSimulation"]["n_neighbors"], dataset=config["datasets"]),
        run_normalization=expand("output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv" ,FC=config["semiSimulation"]["FC"], n_neighbors=config["semiSimulation"]["n_neighbors"], dataset=config["datasets"]),
        run_methods=expand("output/{dataset}/{method}/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.tsv" ,FC=config["semiSimulation"]["FC"], n_neighbors=config["semiSimulation"]["n_neighbors"], dataset=config["datasets"], method=config["methods"], l_param_index=config["l_param_index"]),
        run_metrics=expand("output/{dataset}/metrics/{method}/f1_score_FC_{FC}_n_neigbors_{n_neighbors}_CT_statistics_l_{l_param_index}.csv",FC=config["semiSimulation"]["FC"], n_neighbors=config["semiSimulation"]["n_neighbors"], dataset=config["datasets"], method=config["methods"], l_param_index=config["l_param_index"])
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
        cellmetadata="data/processed/{dataset}/cellmetadata4_{dataset}.json"
    output:
        inflated_counts="output/{dataset}_semiSimulation_NB/inflated_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv",
        simulated_cellmetadata="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}.json",
        simulated_interactions="output/{dataset}_semiSimulation_NB/simulated_interactions_FC_{FC}_n_neigbors_{n_neighbors}.RDS",
        FC_after_simulation="output/{dataset}_semiSimulation_NB/realFC_aftersimulation_{FC}_n_neigbors_{n_neighbors}.RDS",
        metadata_cpdbv5="data/cpdbv5_extrafiles/{dataset}/metadata_cpdbv5_FC_{FC}_n_neigbors_{n_neighbors}.tsv",
        plot_neighbors="data/processed/{dataset}/plot_neighbors_FC_{FC}_n_neigbors_{n_neighbors}.pdf"
    params:
        nLR_per_CTCTcomb=config["nLR_per_CTCTcomb"],
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
        inflated_counts="output/{dataset}_semiSimulation_NB/inflated_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv"
    output:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv"
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
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv",
        microenvironment="data/cpdbv5_extrafiles/microenvironment.tsv",
        metadata="data/cpdbv5_extrafiles/{dataset}/metadata_cpdbv5_FC_{FC}_n_neigbors_{n_neighbors}.tsv",
        cpdb_database="data/cpdbv5_extrafiles/cellphonedb.zip",
    output:
        significant_interactions="output/{dataset}/cellphonedb/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}.tsv"
    container:
        "sing_container/cellphonedbv5.sif"
    script:
        "scripts/method_cellphonedb.py"
        
rule run_method_lianaP_morans:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}.json"
    output:
        significant_interactions="output/{dataset}/lianaP_morans/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.tsv",
        plot_neighbors="output/{dataset}/lianaP_morans/plot_neighbors_lianaP_morans_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.png"
    container:
        "sing_container/lianaPlus.sif"
    params:
      l_index = "{l_param_index}",
      dataset = "{dataset}"
    script:
        "scripts/method_lianaP_morans.py"
        
rule run_method_lianaP_lee:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}.json"
    output:
        significant_interactions="output/{dataset}/lianaP_lee/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.tsv",
        plot_neighbors="output/{dataset}/lianaP_lee/plot_neighbors_lianaP_lee_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.png"
    container:
        "sing_container/lianaPlus.sif"
    params:
      l_index = "{l_param_index}",
      dataset = "{dataset}"
    script:
        "scripts/method_lianaP_lee.py"
        
rule run_method_spatialdm:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}.json"
    output:
        significant_interactions="output/{dataset}/spatialdm/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.tsv",
        plot_neighbors="output/{dataset}/spatialdm/plot_neighbors_spatialdm_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.png"
    container:
        "sing_container/spatialdm.sif"
    params:
      l_index = "{l_param_index}",
      dataset = "{dataset}"
    script:
        "scripts/method_spatialdm.py"
        
        
rule run_method_mistyR:
    threads: 1
    resources:
        mem_mb=5000
    input:
        normalized_counts="output/{dataset}_semiSimulation_NB/inflated_normalized_counts_FC_{FC}_n_neigbors_{n_neighbors}.tsv",
        cellmetadata_post_simulation="output/{dataset}_semiSimulation_NB/simulated_cellmetadata_{FC}_n_neigbors_{n_neighbors}.json"
    output:
        significant_interactions="output/{dataset}/mistyR/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.tsv",
        results_folder=directory("output/{dataset}/mistyR/folder_neighbors_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}")
    params:
    	LR_database=config["LR_database"]
    container:
        "sing_container/mistyR.sif"
    params:
      l_index = "{l_param_index}",
      dataset = "{dataset}"
    script:
        "scripts/method_mistyR.R"
################## METRICS
# f1score
rule run_metric_f1score:
    threads: 1
    resources:
        mem_mb=1000
    input:
        significant_interactions="output/{dataset}/{method}/significant_interactions_FC_{FC}_n_neigbors_{n_neighbors}_l_{l_param_index}.tsv",
        simulated_interactions="output/{dataset}_semiSimulation_NB/simulated_interactions_FC_{FC}_n_neigbors_{n_neighbors}.RDS",
    output:
        CT_statistics="output/{dataset}/metrics/{method}/f1_score_FC_{FC}_n_neigbors_{n_neighbors}_CT_statistics_l_{l_param_index}.csv"
    container:
        "sing_container/liana_edgeR.sif"
    script:
        "scripts/metric_f1score.R"
