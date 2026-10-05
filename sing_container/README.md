# Singularity / Apptainer containers

This folder holds the def files for the containers used by the CCC benchmark pipeline. You can build the images locally from these recipes.

## Contents

```
sing_container/
├── README.md
├── get_def_files.sh      # maintainer script: extracts .def recipes from existing .sif images
└── defs/
    └── <name>.def        # one definition file per container
```

## Requirements

- [Apptainer](https://apptainer.org/docs/admin/main/installation.html) ≥ 1.1 (not tested) **or** [SingularityCE] ≥ 4.0.0. 

## Building the containers

Build a single container:

```bash
cd sing_container
apptainer build --fakeroot <name>.sif defs/<name>.def
# or: sudo apptainer build <name>.sif defs/<name>.def
```

Build all containers in one go:

```bash
cd sing_container
./get_def_files.sh
```

The resulting `.sif` files are written to `sing_container/`, which is where the Snakemake rules expect them.

> **Note:** Builds pull the latest package versions available at build time unless versions are pinned in the `.def` file, so rebuilt images may differ slightly from the ones used in the manuscript.
# Will try soon to pin the package versions for edgeR.sif

## Testing a container

```bash
apptainer exec <name>.sif R --version        # for R-based containers
apptainer exec <name>.sif python --version   # for Python-based containers
```
