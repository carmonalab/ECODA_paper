# Issue running >12h simulations on private-kruse-gpu

- Source: https://hpc-community.unige.ch/t/4291

- Created: 2026-05-04T15:36:29.066Z

- Tags: baobab, slurm

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Alberto.Dinelli (2026-05-04T15:36:29.150Z)

## Primary informations
Username: dinelli
Cluster: baobab
Subject: slurm|private-kruse-gpu
## Description
I have issues running jobs on pivate-kruse-gpu for more than 12 hours.
The bash file to run my jobs can be found in:
```
scratch/HighResolution-Lattice/HighResolution.sh
```
This is its content:
```
#!/bin/env bash
#SBATCH --array=1-7
#SBATCH --partition=private-kruse-gpu
#SBATCH --time=0-72:00:00
#SBATCH --output=%J.out
#SBATCH --mem=3000  
#SBATCH --gpus=1
#SBATCH --constraint=nvidia_a100-pcie-40gb|nvidia_a100_80gb_pcie|nvidia_h100_nvl

module load Julia

cd /srv/beegfs/scratch/users/d/dinelli/HighResolution-Lattice/
srun julia --optimize=3 /home/users/d/dinelli/Code/HighResolution-Lattice/NematoPolar-FFT-Main.jl
```
Note that this is a script I have been using with no issue until last Thursday. Now, when I do sbatch HighResolution.sh, my jobs remain pending with the following justification:
```
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
7884341_[1-7] private-kruse-gpu HighReso  dinelli PD       0:00      1 (ReqNodeNotAvail, Reserved for maintenance)
```
However, when running 12h-simulations:
```
#SBATCH --time=0-12:00:00
```
my jobs run fine. I am surprised since we should be allowed to run jobs for more than 12h on the private partition.
Could you help me understand what is going wrong?
Thanks for the help in advance!
Alberto


## Post 2 by @Adrien.Albert (2026-05-05T08:58:37.892Z)

Dear @Alberto.Dinelli[@Alberto.Dinelli](https://hpc-community.unige.ch/u/alberto.dinelli)
This behavior is expected.
A maintenance on Baobab has been scheduled, and the compute nodes are marked as reserved for maintenance to ensure that no jobs are running when the maintenance starts. As a result, SLURM prevents new jobs from starting if their requested wall time would extend beyond the beginning of the maintenance window.
In your case, the job requests 72 hours, but the remaining time before the maintenance is shorter than that. Therefore, the job stays in PENDING state and cannot start.
You can find the official maintenance announcement here:
Baobab scheduled maintenance: 05-07 May 2026[Baobab scheduled maintenance: 05-07 May 2026](https://hpc-community.unige.ch/t/baobab-scheduled-maintenance-05-07-may-2026/4281) HPC Announce[HPC Announce](https://hpc-community.unige.ch/c/hpc-announce/6)
> Dear users, 
We will be performing software and hardware maintenance on the Baobab HPC cluster from 05 to 07 May 2026, starting at 08:00 +0100. 
You’ll receive an email when the maintenance is finished. The cluster will be completely unavailable during this time, with no access whatsoever (not even to retrieve files). 
We remember you that you our two other HPC clusters are available during the maintenance! 
When submitting a job, make sure that the expected wall time (duration) doesn’t overlap …
Best Regards
