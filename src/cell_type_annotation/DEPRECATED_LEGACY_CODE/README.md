# analysis_JiGL_2026_PBMC_atlas

Cell type annotation of \>1000 PBMC samples with HiTME on an HPC cluster.

------------------------------------------------------------------------

## Workflow Overview

``` text
[ Monolithic h5ad Files on NAS ]
               │
               ▼ (1_prepare_chunks.sh -> 1_prepare_chunks.r)
  [ chunk_1.txt ]  [ chunk_2.txt ]  ...  [ chunk_N.txt ]  (5 samples each)
               │               │
               ▼ (SLURM Array) ▼ (SLURM Array)
         [ Core Job 1 ]  [ Core Job 2 ]  ...              (2.1_process_chunk.sh)
               │               │
               ▼               ▼
     [ sample_1.rds ] ... [ sample_6.rds ] ...            (NAS Output Storage)
```

## Usage Guide

Execute the pipeline steps strictly in the following order on the HPC cluster:


### Download libraries to the repo (e.g. in your home directory)

``` bash
git clone git@github.com:carmonalab/analysis_JiGL_2026_PBMC_atlas.git
cd analysis_JiGL_2026_PBMC_atlas
# If you need to update:
# git fetch origin && git reset --hard origin/master
# Allow execution of script files
chmod +x 1_prepare_chunks.sh 2_submit_hpc_array.sh 2.1_run_worker.sh
```

#### Install pixi to create environment (getting correct R version)

```
cd ~
curl -fsSL https://pixi.sh/install.sh | bash
source ~/.bashrc

# Navigate to your project folder
cd ~/analysis_JiGL_2026_PBMC_atlas

# Initialize a pixi project environment
pixi init

# Add the channels required for data science and single-cell tooling
pixi project channel add conda-forge
pixi project channel add bioconda

# Add R 4.5.2 and system compilation requirements
pixi add r-base=4.5.2 r-jsonlite r-reticulate r-renv r-rcpp compilers zlib libxml2 libpng nlopt cmake pkg-config r-xml2 r-nloptr r-igraph r-leidenbase

# xml2, igraph and some other packages are notoriously difficult to install, so use the pre-compiled ones from pixi
cat << 'EOF' > .Rprofile
# 1. Allow renv to initialize and set up its base sandboxing
if (file.exists("renv/activate.R")) {
    source("renv/activate.R")
}

local({
    # 2. Dynamically resolve the local Pixi library path
    pixi_lib <- file.path(getwd(), ".pixi/envs/default/lib/R/library")
    
    if (dir.exists(pixi_lib)) {
      # Force renv to treat Pixi as a trusted external library link
      options(renv.config.external.libraries = pixi_lib)
      
      # Append the Pixi binary path directly to the R runtime search mapping
      .libPaths(unique(c(.libPaths(), pixi_lib)))
    }
    
    # 3. Dynamically resolve your project's linux-rocky library path
    # R.version$platform will yield 'x86_64-pc-linux-gnu' or similar
    r_mm <- paste0(R.version$major, ".", sub("\\..*$", "", R.version$minor))
    rocky_lib <- file.path(getwd(), "renv/library/linux-rocky-9.6", paste0("R-", r_mm), R.version$platform)
    
    # Fallback just in case the layout skips the platform subfolder
    rocky_lib_fallback <- file.path(getwd(), "renv/library/linux-rocky-9.6")
    
    # 4. Collect all valid targets
    targets <- c(rocky_lib, rocky_lib_fallback, pixi_lib)
    valid_targets <- targets[dir.exists(targets)]
    
    if (length(valid_targets) > 0) {
        # Force renv to treat these directories as trusted external spaces
        options(renv.config.external.libraries = valid_targets)
        
        # Inject them cleanly into the front of the runtime search map
        .libPaths(unique(c(valid_targets, .libPaths())))
    }
})
EOF

# Grab an interactive compute slot to safely compile your packages
srun --partition=shared-cpu --time=01:00:00 --ntasks=1 --cpus-per-task=2 --mem=8G --pty bash

# Run the restore inside your new R 4.5.2 isolation layer
pixi run R -e "renv::install('carmonalab/SignatuR@0aded1db')"
pixi run R -e "renv::restore(repos = c(CRAN = 'https://cloud.r-project.org'))"
pixi run R -e "library('reticulate'); reticulate::install_miniconda(force=TRUE); py_require('anndata')"

# Leave the worker node
exit
```

#### Check and update the config.env file

The config.env file contains all the folder paths and your email for logging.
Adjust accordingly.


### 💡 Crucial Note on Testing & Debugging

Before executing the entire dataset, always perform a quick "micro-test" to ensure environment variables and reference paths are configured correctly:
This runs a small test set of one sample per chunk file and 

``` bash
./1_prepare_chunks.sh test
./2_submit_hpc_array.sh
```

### 0. Download dataset (optional)

Downloads raw single-cell matrices and metadata tables from Zenodo if they are not already present on the storage volume.

Bash

``` bash
Rscript 0_download_dataset.r
```

### 1. Define chunks to be processed on HPC

Bash

``` bash
Rscript 1_prepare_chunks.r
```

-   Input: Raw .h5ad files.

-   Process: Streams target objects in Python backed mode (protecting local system memory) and groups sample arrays into manageable batch files. It also initializes your centralized project paths via get_pipeline_config() and moves files from NAS to HOME on HPC (NAS not accessible from worker nodes).

-   Output: Clean map configurations saved to chunks/chunk_x.txt.

### 2. Submit the HPC Job Array

Bash

``` bash
bash 2_submit_hpc_array.sh
```

-   Process: This launcher automatically reads your configuration file, counts the number of available chunks, and handles your sbatch scheduling metrics smoothly. Each sub-task processes its given samples sequentially on an isolated single CPU core using `2.1_process_chunk.sh`. Also copies back data from HOME to NAS.

-   Output: Sample-level analytical data saved back to your NAS directory structure:

    -   Individual Seurat objects: `output/samples/<sample_id>.rds`

    -   Summary composition lists: `output/ecoda/<sample_id>.rds`

Optional (if ./2_submit_hpc_array.sh failed), copy files manually:

``` bash
rsync -rltvkz --no-p --no-g --progress "/home/users/h/halterc/Jimenez-Gracia_2026_41526507/output/" "/srv/smednas515.unige.ch/carmona_smb/DataCollections/Standardized_SingleCell_Datasets/Jimenez-Gracia_2026_41526507/output/"
```

### 3. Summarize ECODA data

Bash

``` bash
Rscript 3_summarize_ecoda.R
```

-   Input: Collected tracking slices from `output/ecoda/*.rds`.

-   Process: Merges all independent metrics together, scales counts across features, filters low-quality cells, and creates a unified compositional summary object.

-   Output: Centralized master matrix file.
