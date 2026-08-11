# [Baobab] Job not starting while node are idle

- Source: https://hpc-community.unige.ch/t/3770

- Created: 2024-12-19T09:41:20.844Z

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ludovic.Dumoulin (2024-12-19T09:41:20.889Z)

## Primary informations
Username: dumoulil
Cluster: baobab
## Description
My job array is blocked at 26, the 27 is waiting while node of my private gpus are idle
```
private-kruse-gpu    up 7-00:00:00      5   idle gpu[004,006,020,030-031]
```
(Now someone is running jobs on the node GPU020. )
My jobs are not starting due to “ReqNodeNotAvail”:
```
squeue --me
14059875_[27-378%4 private-k C2C_arra dumoulil PD       0:00      1 (ReqNodeNotAvail, UnavailableNodes:gpu[002,013-014,016,027,029,034,036,042-043])
```
while my private nodes are not on the list…
I don’t understand what is happening,
Thank you for your help,
Best


## Post 2 by @Ludovic.Dumoulin (2024-12-19T09:44:51.349Z)

If I want to submit the same bash again, but only on private gpus, it is not working:
> sbatch C2C_array.sh
> sbatch: error: Batch job submission failed: Requested node configuration is not available
My bash:
```
#!/bin/env bash
#SBATCH --array=1-378%40
#SBATCH --partition=private-kruse-gpu
#SBATCH --time=0-12:00:00
#SBATCH --output=%J.out
#SBATCH --mem=3000  
#SBATCH --gpus=ampere:1 
#SBATCH --constraint=DOUBLE_PRECISION_GPU

module load Julia

cd /srv/beegfs/scratch/users/d/dumoulil/Data/P-series/2defects_det/
srun julia --optimize=3 /home/users/d/dumoulil/Code/FFT_P_2def_det/2D.jl
```


## Post 3 by @Yann.Sagon (2024-12-19T13:03:37.605Z)

Dear @Ludovic.Dumoulin[@Ludovic.Dumoulin](https://hpc-community.unige.ch/u/ludovic.dumoulin)
Please check the announcement I just posted, you’ll find the reason. Important: new GPU types naming in Baobab[Important: new GPU types naming in Baobab](https://hpc-community.unige.ch/t/important-new-gpu-types-naming-in-baobab/3772)
