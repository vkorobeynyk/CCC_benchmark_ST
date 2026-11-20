# CCC_benchmark_ST

A reproducible benchmark of **spatial cell-cell communication (CCC) methods** on spatial transcriptomics data.

Real spatial datasets are **semi-simulated**: the native expression of a known ligand-receptor (LR) pair is replaced by simulated signal in chosen Sender and Receiver cells. Each method is then scored on whether it recovers that pair. Each method is also run on the original data without simulated signal (**baseline**), so interactions caused by the simulation can be separated from those already present in the data. The workflow is built with [Snakemake](https://snakemake.readthedocs.io/), runs on a SLURM cluster, and runs every step inside a Singularity/Apptainer container.

This repository is the spatial part of a two-part benchmark. The single-cell part lives in the companion `CCC_benchmark` repository.

The data folder used in the benchmark, along with the output folder can be found in ................

---

## Contact

Vladyslav Korobeynyk, HIFO / DMLS, University of Zurich (Jessberger lab / Mark D. Robinson lab)
Questions and bug reports: please open a [GitHub issue](https://github.com/vkorobeynyk/CCC_benchmark_ST/issues).

## Generative AI statement
Generative AI was used throughout this entire benchmark to make code nicer to read and more efficient. The entire logic of the benchmark was created by myself and I assume responsability of the content within this repo.

## Citation

<!-- TODO: add manuscript reference / preprint link -->

## License

This project is released under the [GNU General Public License v3.0](LICENSE).


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

Seurat Wilcoxon has no spatial information. It is the baseline for measuring how much the spatial methods gain from using space. `Table2_spatial_Cell_Communication_Summary.pdf` (made by `generate_table.R`) summarises how each method scores interactions and tests significance.

## Datasets

Set in `config.yaml` → `datasets`:

| ID | Technology | Tissue | Source |
|---|---|---|---|
| `Visium_HD_HPC` | 10x Visium HD (binned to cells with Bin2cell) | Mouse brain / hippocampus | [10x Genomics](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-mouse-brain-fresh-frozen) |
| `MERFISH_mColon` | MERFISH | Mouse colon | [Dryad](https://datadryad.org/dataset/doi:10.5061/dryad.rjdfn2zh3) (Cadinu et al., 2024, *Cell*) |
| `CosMx_HFC` | NanoString CosMx | Human frontal cortex (FFPE) | [NanoString](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-frontal-cortex-ffpe-dataset/) |

Raw data is **not** in the repository (`data/` is git-ignored). See [Input data](#input-data) for the expected layout. Download notes are in `data/download_data.R`.

### Dataset preprocessing

Before entering the pipeline, each dataset was converted into a cell × gene count matrix (`counts_<dataset>.tsv`) and a per-cell metadata table (`metadata_<dataset>.tsv` with `Cell_ID`, `Celltype`, `x`, `y`).

- **`Visium_HD_HPC`**: The 2 µm bins were assigned to cells with [Bin2cell](https://github.com/Teichlab/bin2cell), and counts were summed per cell to give a gene × cell matrix.
- **`CosMx_HFC`**: Coordinates were converted from millimetres to micrometres. Some cell types were removed, reducing the original ~190k cells to ~32k. That set was then randomly subsampled to 1/4 (~8k cells) to lower computational cost.

---

## Benchmark design

### Two semi-simulation strategies

Both strategies run in a single Snakemake call and write to separate output trees (`output/<strategy>/…`, `data/processed/<strategy>/…`).

| Strategy | How Sender/Receiver cells are chosen | Question it asks |
|---|---|---|
| `spatialScattering` | Randomly across the whole tissue, so Sender-Receiver distances cover the dataset's full range | Can methods find a signal that is spread out in space? |
| `spatialColocalization` | From a hand-picked group of cells that sit together in one region (chosen with `shiny_script_forSelectingCells.R`), split at random into Sender and Receiver halves | Can methods find a signal that is spatially concentrated? |

### Semi-simulation

For each dataset, the pipeline:

1. Estimates per-gene mean and dispersion with **edgeR**, using real per-cell offsets.
2. Takes **one LR pair** from `data/LR_database.tsv` (`indexLR_toSample`). Pairs are simulated one at a time so the already sparse spatial data is not distorted.
3. Replaces the native expression of the pair's genes: the ligand is first set to 0 in all Sender cells and the receptor in all Receiver cells, then re-expressed in the selected cells only. New counts are drawn from a negative binomial, `1 + rnbinom(mu = gene_rate · exp(cell_offset), size = 1/dispersion)`. For multi-subunit complexes, every subunit gets the signal. A fixed seed per parameter combination (`PCE_Sender`, `PCE_Receiver`, `indexLR_toSample`) makes the cell selection and sampling reproducible.
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

`radius_param` in `config.yaml` sets the neighbourhood size (Euclidean radius or kernel bandwidth) for **each method and dataset**, in that dataset's coordinate units.

### Baseline runs (no simulated signal)

Every method is also run on the original processed data, without semi-simulation (`scripts/processing_baseline.R` and the `*_baseline` rules in the `Snakefile`). Methods infer interactions from the fixed `Celltype` labels (Sender/Receiver), so the baseline does not depend on `PCE_Sender`, `PCE_Receiver` or `indexLR_toSample`. It runs once per **strategy × dataset × radius_param_index × method** (144 runs with the default config).

- The method rules are reused with Snakemake's `use rule ... with:`, so scripts and containers are identical to the semi-simulated runs.
- `Celltype_updated` is only used by the methods for diagnostic columns (e.g. how many signal-added cells fall within the radius). In the baseline, **all** Sender cells are labelled `Sender_signalAdded` and all Receiver cells `Receiver_signalAdded`, so these diagnostics run without errors. Baseline diagnostic columns are not used downstream.

### Scoring

`scripts/metric_f1score_rankingLRgenes.R` gathers every method's `significant_interactions_*.tsv` for one strategy. It treats the simulated Sender → Receiver LR pair as ground truth and compares each run with the baseline of the same strategy, dataset, method and radius:

- **TP / FN**: whether the simulated pair is called significant.
- Every other significant interaction falls into exactly one group:
  - **pre-existing**: also significant in the baseline. Present in the data regardless of the simulation; not counted as FP.
  - **related**: not significant in the baseline, but shares a ligand, receptor or complex subunit with the simulated pair; not counted as FP.
  - **FP**: not significant in the baseline and shares no gene with the simulated pair, i.e. a spurious hit caused by the simulation.
- **Precision, recall, F1** are computed from TP, FN and this FP, together with the **rank** of the simulated pair.
- Additional columns per run:
  - `n_additional`: all significant non-simulated interactions (`n_pre_existing + n_related + FP`); the FP definition before baseline runs were added
  - `n_pre_existing`, `n_related`: sizes of the groups above
  - `n_lost`: significant in the baseline but not after simulation, split into `n_lost_simulated`, `n_lost_related` and `n_lost_unrelated`. Because the simulation removes the native expression of the simulated genes, simulated and related losses are expected; unrelated losses would indicate a global effect (e.g. normalization or cross-pair scaling).
  - `n_baseline_significant`: total significant interactions in the baseline run

Results are saved to `output/<strategy>/final_scores.RDS`, and `scripts/visualization.R` turns them into the summary figures. Each `<dataset>_recall_precision_plots.pdf` ends with two baseline-comparison plots per radius scenario: the fate of baseline-significant interactions, and the composition of non-simulated significant interactions (pre-existing / related / FP).

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
2. run_semiSimulation_inflateCounts   plant the LR signal  (processing_semisimulation_nbsampling.R)
        │
        ▼
3. run_normalization                  scater::logNormCounts  (processing_normalization.R)
        │
        ▼
4. run_method_*                       8 CCC methods on the same normalised input
        │
        │  4b. baseline (no simulated signal), per strategy × dataset × radius:
        │      run_baseline_prepare           original counts + metadata  (processing_baseline.R)
        │      run_baseline_normalization     scater::logNormCounts
        │      run_method_*_baseline          same 8 methods  → <dataset>/<method>/baseline/
        │
        ▼
5. run_metric_f1score_rankingLRgenes_all   signal runs vs. baseline → output/<strategy>/final_scores.RDS
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
├── fix_indentation.sh                # normalises indentation in Snakefile/config before a run
├── generate_table.R                  # builds Table 2 (method summary)
├── shiny_script_forSelectingCells.R  # interactive cell picker for spatialColocalization
├── merge_svg_to_pdf.sh               # figure post-processing helper
├── scripts/
│   ├── processing_dataset_spatialScattering.R
│   ├── processing_dataset_spatialColocalization.R
│   ├── processing_semisimulation_nbsampling.R
│   ├── processing_baseline.R
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

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 6 (the baseline rules use `use rule ... with:`) with a SLURM profile (`--profile slurm`)
- [Apptainer](https://apptainer.org/) or [SingularityCE](https://sylabs.io/singularity/) ≥ 4.0
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
│   ├── metadata_<dataset>.tsv               # per-cell metadata incl. spatial coordinates
│   └── selected_cells.csv                   # spatialColocalization only (Cell_ID column);
│                                            #   Visium_HD_HPC uses DG_cells.csv
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
snakemake --use-singularity --profile slurm output/spatialScattering/visualize_results.done
```

## Outputs

```
output/<strategy>/
├── <dataset>_semiSimulation_NB/     # inflated + normalised counts, simulated metadata, ground truth
├── <dataset>_baseline/              # original + normalised counts and metadata, no simulated signal
├── <dataset>/<method>/              # significant_interactions_*.tsv and neighbourhood plots
├── <dataset>/<method>/baseline/     # significant_interactions_baseline_l_<radius>.tsv (baseline runs)
├── final_scores.RDS                 # precision / recall / F1 / rank and baseline-comparison counts per combination
└── …                                # summary figures from visualization.R
data/processed/<strategy>/<dataset>/ # processed counts, edgeR estimates, QC and distance plots
data/cpdbv5_extrafiles/<strategy>/<dataset>/baseline/  # CellPhoneDB metadata for the baseline
```

---

## Notes

- **stLearn** runs through `shell:` with an argparse CLI, not `script:`. Snakemake's script mode clashes with stLearn's Python 3.10 environment.
- **CellChat** sets `scale.distance` from the data (`1.5 / min_pairwise_distance`) so it works across datasets with different coordinate scales. This changes only the raw probability values, not significance or ranks.
- Container builds install the latest package versions unless versions are pinned in the `.def` files, so rebuilt images may differ slightly from the ones used for the manuscript.
**TODO: add versions to .def**
