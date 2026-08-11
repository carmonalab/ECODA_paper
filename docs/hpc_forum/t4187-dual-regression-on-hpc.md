# Dual regression on HPC

- Source: https://hpc-community.unige.ch/t/4187

- Created: 2026-01-07T15:12:06.810Z

- Tags: baobab

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Faustine.Lemarechal (2026-01-07T15:12:06.860Z)

Hello, i am trying to do a dual regression on HPC, but it’s not working. Each time it seems to work but I don’t have any output.
Here is my code:
```
module load GCC/8.3.0

module load OpenMPI/3.1.4

module load FSL/6.0.3-Python-3.7.4

\# The melodic output

MELODIC_IC="/home/users/l/lemarec4/data_ICA/GroupICA_T0_T1.ica/melodic_IC.nii.gz"

\# Output directory for dual regression results

RESULTS_DIR="/home/users/l/lemarec4/Results_dualreg"

\# Must be the SAME list (same order) used for MELODIC

INPUT_LIST="/home/users/l/lemarec4/data_ICA/inputs_4_ica.txt"

\# Perform thresholded dual regression to obtain unbiased timeseries for connectomics analyses (e.g., with FSLnets)

THR_FLAG="--thr"

mkdir -p "${RESULTS_DIR}"

\# ---------------------------

\# 3) Run Dual Regression (Stages 1 & 2 only)

\# ---------------------------

\# Arguments: <groupICs> <des_norm> <design.mat> <design.con> <n_perm> <outdir> <subjects...>

\# Here: "0 -1 0" means: no group-level stage 3 here (stages 1 & 2 only)

echo "Running dual_regression on inputs from: ${INPUT_LIST}"

srun dual_regression   "${MELODIC_IC}"   0 -1 0 ${THR_FLAG}   "${RESULTS_DIR}"   \`cat "${INPUT_LIST}"\`

echo "Done. Main outputs:"

echo " - ${RESULTS_DIR}/dr_stage1/ (time courses)"

echo " - ${RESULTS_DIR}/dr_stage2_icXXXX.nii.gz (subject spatial maps per IC)"
```
Can you help me ?
Thank you so much !


## Post 2 by @Yann.Sagon (2026-01-08T16:24:20.978Z)

Dear @Faustine.Lemarechal[@Faustine.Lemarechal](https://hpc-community.unige.ch/u/faustine.lemarechal)
Faustine.Lemarechal:
```
module load GCC/8.3.0

module load OpenMPI/3.1.4

module load FSL/6.0.3-Python-3.7.4
```
Please try with the latest version, this one is very old.
```
ml GCC/10.3.0  OpenMPI/4.1.1 FSL/6.0.5.1
```
Faustine.Lemarechal:
> Each time it seems to work but I don’t have any output.
Which command do you use to start your job?


## Post 3 by @Faustine.Lemarechal (2026-01-17T10:36:36.745Z)

Hello, thank you so much for your help. But it didn’t work even with your advices. It seems that BAOBAB does not include fsl_sub which the script includes. There is my script:
```
#!/bin/bash

#SBATCH --time=10:00:00               # NOTE: likely longer than generally needed 

#SBATCH --ntasks 1

#SBATCH --cpus-per-task=16

#SBATCH --partition=shared-cpu

\# Outputs ----------------------------------

#SBATCH --output log/%x-%A-%a.out

#SBATCH --error log/%x-%A-%a.err

#SBATCH --mail-user=faustine.lemarech@hes-so.ch

#SBATCH --mail-type=ALL

\# ------------------------------------------

\# ---------------------------

\# USER PARAMETERS TO EDIT

\# ---------------------------

MELODIC_IC="/home/users/l/lemarec4/data_ICA/GroupICA_T0_T1.ica/melodic_IC.nii.gz"

RESULTS_DIR="/home/users/l/lemarec4/Results_dualreg"

INPUT_LIST="/home/users/l/lemarec4/data_ICA/inputs_4_ica.txt"

THR_FLAG="--thr"

\# Crée le dossier de sortie si nécessaire

mkdir -p "${RESULTS_DIR}"

\# ---------------------------

\# Load FSL (Baobab modules)

\# ---------------------------

ml purge

ml GCC/10.3.0 OpenMPI/4.1.1 FSL/6.0.5.1

\# Permet à FSL de trouver les packages Python (nibabel, etc.)

export PYTHONPATH=/home/users/l/lemarec4/.local/lib/python3.9/site-packages:$PYTHONPATH

export FSLPARALLEL=0

\# ---------------------------

\# Run Dual Regression (Stages 1 & 2 only)

\# ---------------------------

echo "Running dual_regression on inputs from: ${INPUT_LIST}"

dual_regression "${MELODIC_IC}" 0 -1 0 ${THR_FLAG} "${RESULTS_DIR}" $(cat "${INPUT_LIST}") --nobet

echo "Done. Main outputs:"

echo " - ${RESULTS_DIR}/dr_stage1/ (time courses)"

echo " - ${RESULTS_DIR}/dr_stage2_icXXXX.nii.gz (subject spatial maps per IC)"
```
Thank you again


## Post 4 by @Yann.Sagon (2026-01-19T09:54:06.364Z)

Faustine.Lemarechal:
> `# Permet à FSL de trouver les packages Python (nibabel, etc.)`
You can load the module `NiBabel/3.2.1` instead.
Faustine.Lemarechal:
> `export FSLPARALLEL=0`
This is for testing? For production later, be careful to specify the number of cores you request.
Faustine.Lemarechal:
> It seems that BAOBAB does not include fsl_sub which the script includes.
I’ll install a newer version, it seems it is now installed by default.


## Post 5 by @Yann.Sagon (2026-01-19T16:43:07.564Z)

New software installed: FSL version 6.0.7.17[New software installed: FSL version 6.0.7.17](https://hpc-community.unige.ch/t/new-software-installed-fsl-version-6-0-7-17/4199) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Dear users, we have installed a new software: FSL 6.0.7.17: 

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  FSL: FSL/6.0.7.17
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------…
