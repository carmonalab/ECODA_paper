see datasets.json for most up-to-date list of datasets used and conditions

# HPC
- bash SLURM submission scripts are run on the login node, spawning worker nodes
- only login node has access to the shared NAS file system
- worker nodes do NOT have access to NAS
- data must be copied to local scratch before processing (done with ./src/bash/copy_data_from_nas_to_hpc_scratch.sh)
- results must be copied back to NAS after processing
- If more information is needed, documentation for the HPC can be found here: https://doc.eresearch.unige.ch/hpc/start