# Paper/repo review and update strategy
This repo is about: ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts
Link to paper: https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full

## Major reviewer points to be addressed
- Extend batch effect analysis
    - More datasets have to be added
        - A requirement for this: Batch effect analysis pipeline has to be implemented in a more structured and standardized way
        - Preprocessing has to be harmonized to fit with the benchmark analysis
        - Additional datasets will be added by me (human) and do not need to be addressed by agent
    - More methods have to be added
        - Methods to be added are defined in the TODO.md
        - A draft strategy for specific pipeline and code implementation are defined in the TODO.md
- Extend benchmark analysis
    - More datasets have to be added
        - Benchmark analysis has mainly to be adapted to be run on the HPC cluster and to be cleaned up with minor adaptions
        - Additional datasets will be added by me (human) and do not need to be addressed by agent
    - More methods have to be added
        - Methods to be added are PILOT-GM-VAE (very similar to PILOT which is already implemented, trivial) and possibly PULSAR (needs to be tested if it can be run at all, as it is a foundation model, requiring substantial hardware only available on the cluster, including GPU. Input for PULSAR is Universal Cell Embedding (UCE)).
            - PILOT-GM-VAE can be added by agent
            - PULSAR needs to be tested for requirements


# datasets.json
This acts as ground truth for the datasets evaluated in this study. See datasets.json for most up-to-date list of datasets used and conditions.
Do not change this file without asking.

# Batch effect analysis dataset info
- Most datasets are monolithic h5ad files with a batch_col, e.g.:
    - Stephenson has batch effect by batch_col "Site" (both, in terms of gene expression (major) and cell type composition (minor, just one or a few monocyte subtypes))
- A "combined PBMC" dataset was created from multiple other available datasets (included for method benchmark analysis and or batch effect analysis):
    - Combined PBMC (Stephenson, GongSharma, Zhu) (see batch_effect_analysis.rmd, see also TODO.md)


# HPC
- bash SLURM submission scripts are run on the login node, spawning worker nodes
- only login node has access to the shared NAS file system
- worker nodes do NOT have access to NAS
- data must be copied to local scratch before processing (done with ./src/bash/copy_data_from_nas_to_hpc_scratch.sh)
- results must be copied back to NAS after processing
- If more information is needed, documentation for the HPC can be found here: https://doc.eresearch.unige.ch/hpc/start