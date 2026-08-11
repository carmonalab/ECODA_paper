# Sbatch running error with public-cpu

- Source: https://hpc-community.unige.ch/t/4015

- Created: 2025-07-25T08:01:00.470Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jade.Awada (2025-07-25T08:01:00.512Z)

Dear HCP Team,
On Baobab (public-cpu), I’m trying to run a MATLAB script with parallel processing via sbatch, but I’m having some issues, and I am wondering if the configuration parameters are valid or not.
When I try to launch the script, I get this message: sbatch: error: Batch job submission failed: Requested node configuration is not available
Please find below the parameters:
```
#!/bin/bash
#SBATCH --job-name=split_half
#SBATCH --output=split_half_%j.out
#SBATCH --error=split_half_%j.err
#SBATCH --partition=public-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # how many cores MATLAB should use
#SBATCH --mem-per-cpu=10G # memory per core
#SBATCH --time=4-00:00:00 # wall‑time hh:mm:ss
#SBATCH --mail-type=BEGIN,END,FAIL

echo “=== JOB START: $(date) ===”
module load MATLAB/2022a

srun matlab -nodisplay -r “run(‘/home/users/a/awada/consensus_clustering_CAPs_deactivated_GVA_T1/e_Script_SplitHalf_ consensus_clustering_PARALLEL_SLURM.m’); exit;”

echo “=== JOB END: $(date) ===”
```
Can you please help or advice on how to best configure the parameters for running this?
many thanks in advance,
Jade


## Post 2 by @Yann.Sagon (2025-07-29T06:41:55.092Z)

Dear @Jade.Awada[@Jade.Awada](https://hpc-community.unige.ch/u/jade.awada) on Baobab cluster, the public-cpu partition is very small yet and the nodes are old: they have 20 cores with 96GB of ram.
```
(baobab)-[root@admin1 ~]$ sinfo -p public-cpu -N -o "%n %c %m"
HOSTNAMES CPUS MEMORY
cpu237 20 96000
cpu238 20 96000
cpu239 20 96000
cpu240 20 96000
cpu241 20 96000
cpu242 20 96000
cpu243 20 96000
cpu244 20 96000
```
You are requesting 16 cores with 10G per core => 160GB which is more than the compute nodes have.
Your options are:
- use less ram
- use shared-cpu partition which is bigger and have more powerfull nodes
- use Bamboo cluster which is newer, public-cpu have 128 cores and 512GB ram.
