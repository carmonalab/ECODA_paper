# [SOLVED] Permission denied when launching a bash script

- Source: https://hpc-community.unige.ch/t/3455

- Created: 2024-05-22T12:34:25.772Z

- Posts: 1

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Esther.SoleMarti (2024-05-22T12:34:25.824Z)

Edit:
there is a typo in the first line of my script, a space between `#!/bin/` and `bash`. Sorry for the inconvenience!
---
Hello all,
I’m having trouble when using a bash script on Baobab. More on it in the description.
## Primary informations
Username: solemart
Cluster: Baobab
## Description
After logging in to Baobab via SSH, I launch a bash script with the command `sbatch SCRIPT_NAME.sh`. I am allocated a job ID but when I check the output file it only has the following error message: `slurmstepd: error: execve(): /var/spool/slurmd/job10182489/slurm_script: Permission denied`.
A search online[search online](https://bugs.schedmd.com/show_bug.cgi?id=2180) about this error message suggests that it may be a problem on my side, but I’m unfamiliar with everything else explained in the linked page and don’t dare to touch anything before asking.
## Steps to Reproduce
- ssh solemart@baobab2.hpc.unige.ch
- sbatch bash_script.sh (script below)
### Script
`#!/bin/ bash
#SBATCH --partition=shared-gpu
#SBATCH --time=04:00:00
#SBATCH --gpus=1
#SBATCH --output=kraken-%j.out
#SBATCH --mem=48GB
#SBATCH --ntasks=12
module purge
module load CUDA/11.8.0 GCCcore/11.2.0 Python/3.9.6
source ~/kraken-env/bin/activate
echo “ManuMcFondue finetuning”
srun ketos train -f alto -d cuda:0 -B 8 -r 0.0003 -u NFC --workers 12 --augment --precision 16 --lag 5 -i RV_train/data/ManuMcFondue.mlmodel --resize union -t RV_train/data/split/train.txt -e RV_train/data/split/eval.txt
echo “Araucania finetuning”
srun ketos train -f alto -d cuda:0 -B 8 -r 0.0003 -u NFC --workers 12 --augment --precision 16 --lag 5 -i RV_train/data/HTR-Araucania_XIX.mlmodel --resize union -t RV_train/data/split/train.txt -e RV_train/data/split/eval.txt`
## Expected Result
I would expect two finetuning jobs to be performed and a set of models to be saved into my profile on Baobab.
## Actual Result
The `.out` file only contains the message: `slurmstepd: error: execve(): /var/spool/slurmd/job10182489/slurm_script: Permission denied`
Thank you for your time & help!
