# CCC_benchmark_ST

A reproducible benchmark of **spatial cell-cell communication (CCC) methods** on spatial transcriptomics data.

Real spatial datasets are **semi-simulated**: a known ligand-receptor (LR) signal is added to selected Sender and Receiver cells. Each method is then scored on whether it recovers that pair. The workflow is built with [Snakemake], runs on a SLURM cluster, and runs every step inside a Singularity container.

This repository is the spatial part of a two-part benchmark. The single-cell part lives in the `CCC_benchmark` repository.

---

## Contact

Vladyslav Korobeynyk, HIFO / DMLS, University of Zurich (Jessberger lab / Mark D. Robinson lab)
Questions and bug reports: please open a [GitHub issue](https://github.com/vkorobeynyk/CCC_benchmark_ST/issues).

## Generative AI statement
Generative AI was used throughout this entire benchmark to make code nicer to read and more efficient. The entire logic of the benchmark was created by myself and I assume responsability of the content within this repo.

## Citation 
<!-- TODO: add manuscript reference / preprint link -->

## Methods benchmarked

| Method | Language | Container | Spatially aware |
|---|---|---|---|
| [CellChat](https://github.com/jinworks/CellChat) | R | `cellchat.sif` | ✅ |
| [NICHES](https://github.com/msraredon/NICHES) | R | `NICHES.sif` | ✅ |
| [LIANA+](https://github.com/saezlab/liana-py) (bivariate Moran's I) | Python | `lianaPlus.sif` | ✅ |
| [MISTy](https://github.com/saezlab/mistyR) (via LIANA+) | Python | `lianaPlus.sif` | ✅ |
| [SpatialDM](https://github.com/StatBiomed/SpatialDM) | Python | `spatialdm.sif` | ✅ |
| [stLearn](https://github.com/BiomedicalMachineLearning/stLearn) | Python | `stlearn.sif` | ✅ |
| [CellPhoneDB v5](https://github.com/ventolab/CellphoneDB) (with microenvironments) | Python | `cellphonedbv5.sif` | partly (microenvironment file) |
| Seurat Wilcoxon (baseline) | R | `liana_edgeR.sif` | ❌ |

Seurat Wilcoxon has no spatial information. It is the baseline for measuring how much the spatial methods gain from using spacial coordinates. `Table2_spatial_Cell_Communication_Summary.pdf` (made by `generate_table.R`) summarises how each method scores interactions and tests significance.

## Datasets

Set in `config.yaml` → `datasets`:

| ID | Technology | Tissue |
|---|---|---|
| `Visium_HD_HPC` | 10x Visium HD (binned to cells with Bin2cell) | Mouse brain / hippocampus |
| `MERFISH_mColon` | MERFISH | Mouse colon |
| `CosMx_HFC` | NanoString CosMx | <!-- TODO: tissue --> |

Raw data is **not** in the repository (`data/` is git-ignored). See [Input data](#input-data) for the expected layout. Download notes are in `data/download_data.R`.

---

## Benchmark design

### Two semi-simulation strategies

Both strategies run in a single Snakemake call and write to separate output trees (`output/<strategy>/…`, `data/processed/<strategy>/…`).

| Strategy | How Sender/Receiver cells are chosen | Question it asks |
|---|---|---|
| `spatialScattering` | Randomly across the whole tissue, so Sender-Receiver distances cover the dataset's full range | Can methods find a signal that is spread out in space? |
| `spatialColocalization` | From a hand-picked group of cells that sit together in one region (chosen with `shiny_script_forSelectingCells.R`), split at random into Sender and Receiver halves | Can methods find a signal that is spatially colocalized? |

### Semi-simulation

For each dataset, the pipeline:

1. Estimates per-gene mean and dispersion with **edgeR**, using real per-cell offsets.
2. Takes **one LR pair** from `data/LR_database.tsv` (`indexLR_toSample`). Pairs are simulated one at a time so the already sparse spatial data is not distorted.
3. Increases ligand counts in some Sender cells and receptor counts in some Receiver cells. New counts are drawn from a negative binomial, `1 + rnbinom(mu = gene_rate · exp(cell_offset), size = 1/dispersion)`. For multi-subunit complexes, every subunit gets the signal.
4. Sets how many cells get the signal with `PCE_Sender` / `PCE_Receiver`:

   ```
   n_cells_inflated = avg_pct_cells_expressing_LR × PCE × n_Sender_or_Receiver_cells
   ```

   > ⚠️ Here `PCE_*` are **multipliers** of the dataset's average expression rate. They are **not** percentages of expressing cells, as in the single-cell benchmark.

### Parameter grid

Each `(strategy, dataset, method)` combination is run over:

| Wildcard | Meaning | Default values |
|---|---|---|
| `PCE_Sender` | Multiplier for the number of Sender cells given the signal | 1, 2, 4, 6, 8, 10 |
| `PCE_Receiver` | Multiplier for the number of Receiver cells given the signal | 1, 2, 4, 6, 8, 10 |
| `indexLR_toSample` | Which LR pair is simulated | 1–10 |
| `radius_param_index` | Which neighbourhood size is used (index into `radius_param`) | 0, 1, 2 |

`radius_param` in `config.yaml` sets the neighbourhood size (Euclidean radius or kernel bandwidth) for **each method and dataset**, in that dataset's coordinate units. The values in `config.yaml` were selected by trial and error (as methods calculate distances differently) to keep the similar amount of spatial neighbors.

### Scoring

`scripts/metric_f1score_rankingLRgenes.R` gets every method's `significant_interactions_*.tsv` for one strategy. It treats the simulated Sender → Receiver LR pair as ground truth and computes **precision, recall, F1 and the rank** of the true pair. Results are saved to `output/<strategy>/final_scores.RDS`, and `scripts/visualization.R` turns them into the summary figures.

---

## Pipeline overview

```
data/<dataset>/{counts,metadata}_<dataset>.tsv
        │
        ▼
1. run_processing                     QC plots, Sender/Receiver selection, neighbour pairing,
                                      edgeR mean/dispersion  (processing_dataset_<strategy>.R)
        │
        ▼
2. run_semiSimulation_inflateCounts   add the LR signal  (processing_semisimulation_nbsampling.R)
        │
        ▼
3. run_normalization                  scater::logNormCounts  (processing_normalization.R)
        │
        ▼
4. run_method_*                       8 CCC methods on the same normalised input
        │
        ▼
5. run_metric_f1score_rankingLRgenes_all   → output/<strategy>/final_scores.RDS
        │
        ▼
6. visualize_results                  → output/<strategy>/ figures
```

## Repository structure

```
CCC_benchmark_ST/
├── Snakefile                         # workflow definition (rules for steps 1–6)
├── config.yaml                       # datasets, strategies, methods, parameter grids
├── sbatch_submit                     # SLURM submission script for the Snakemake driver
├── fix_indentation.sh                # normalises indentation in Snakefile/config before a run (probably just local problem)
├── generate_table.R                  # builds Table 2 (method summary)
├── shiny_script_forSelectingCells.R  # interactive cell picker for spatialColocalization
├── scripts/
│   ├── processing_dataset_spatialScattering.R
│   ├── processing_dataset_spatialColocalization.R
│   ├── processing_semisimulation_nbsampling.R
│   ├── processing_normalization.R
│   ├── general_spatial_QCpipeline_R_function.R
│   ├── helper_functions.R
│   ├── method_cellchat.R
│   ├── method_NICHES.R
│   ├── method_seurat_wilcoxon.R
│   ├── method_cellphonedbv5.py
│   ├── method_lianaP_morans.py
│   ├── method_mistyR.py
│   ├── method_spatialdm.py
│   ├── method_stlearn.py
│   ├── metric_f1score_rankingLRgenes.R
│   └── visualization.R
├── sing_container/
│   ├── README.md                     # how to build the containers
│   ├── get_def_files.sh
│   └── defs/*.def                    # Apptainer/Singularity recipes
├── data/                             # (git-ignored) inputs, see below
└── output/                           # (git-ignored) results
```

---

## Getting started

### Requirements

- [Snakemake](https://snakemake.readthedocs.io/) with a SLURM profile (`--profile slurm`)
- [Apptainer](https://apptainer.org/) (not tested) or [SingularityCE](https://sylabs.io/singularity/) ≥ 4.0
- A SLURM cluster (the pipeline can also run locally if you drop `--profile slurm`, but the full grid is large)

### 1. Build the containers

Recipes are in `sing_container/defs/`. See [`sing_container/README.md`](sing_container/README.md) for full details.

```bash
cd sing_container
apptainer build --fakeroot cellchat.sif defs/cellchat.def
# …repeat for: cellphonedbv5, giotto, liana_edgeR, lianaPlus, mistyR, NICHES, spatialdm, stlearn
```

The Snakemake rules look for the `.sif` images in `sing_container/`.

### 2. Input data

Each dataset in `config.yaml` needs:

```
data/
├── LR_database.tsv                          # curated LR database (see below)
├── <dataset>/
│   ├── counts_<dataset>.tsv                 # genes × cells raw count matrix
│   ├── metadata_<dataset>.tsv               # per-cell metadata incl. spatial coordinates, Sender and Receiver celltypes
│   └── selected_cells.csv                   # spatialColocalization only (Cell_ID column)
└── cpdbv5_extrafiles/
    ├── microenvironment.tsv                 # CellPhoneDB microenvironment definition
    └── cellphonedb_<version>.zip            # CellPhoneDB v5 database
```

**LR database.** `data/LR_database.tsv` comes from CellPhoneDB in LIANA format. Complex subunits are joined with `_` (e.g. `L_R1_R2`). Self-interactions and genes with a dash in their name were removed, and only pairs also in `CellChatDB.human` were kept. The rows were then shuffled, because the simulation picks LR pairs by row index. The full curation code is in `data/InfoAbout_LR_database.txt`.

> Visium HD has only a few multi-subunit LR pairs, so keep `indexLR_toSample` below 16 for that dataset.

### 3. Configure

Edit `config.yaml` to choose datasets, strategies, methods and parameter grids. In `sbatch_submit`, change the `--bind` path to your own working directory.

### 4. Run

```bash
# dry run
snakemake -n --use-singularity

# full run on SLURM
sbatch sbatch_submit
```

To run one strategy only:

```bash
snakemake --use-singularity output/spatialScattering/visualize_results.done
```

## Outputs

```
output/<strategy>/
├── <dataset>_semiSimulation_NB/     # inflated + normalised counts, simulated metadata, ground truth
├── <dataset>/<method>/              # significant_interactions_*.tsv and neighbourhood plots
├── final_scores.RDS                 # precision / recall / F1 / rank per combination
└── …                                # summary figures from visualization.R
data/processed/<strategy>/<dataset>/ # processed counts, edgeR estimates, QC and distance plots
```

---

## Notes

- **stLearn** runs through `shell:` with an argparse CLI, not `script:`. Snakemake's script mode clashes with stLearn's Python 3.10 environment.
- **CellChat** sets `scale.distance` from the data (`1.5 / min_pairwise_distance`) so it works across datasets with different coordinate scales. This changes only the raw probability values, not significance or ranks.
- Container builds install the latest package versions unless versions are pinned in the `.def` files, so rebuilt images may differ slightly from the ones used for the manuscript.
**To do: Have to add the version for the last point**
