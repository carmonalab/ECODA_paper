# Random: Bus error

- Source: https://hpc-community.unige.ch/t/3249

- Created: 2024-01-12T17:02:32.186Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @parzival.nussbaum (2024-01-12T17:02:32.247Z)

Username: nussbaup
Cluster: Yggdrasil
Launching a job array yields in failing jobs due to a “Bus error”.
I tried to launch this batch file multiple times. The issue seems neither not to be limited to a single node, nor is limited to a specific task in the array.
Is this issue linked to the current /scratch issue on Yggdrasil?
## Steps to Reproduce
Launch the following batch file.
```
#!/bin/bash
#SBATCH --job-name infere_validation    # this is a parameter to help you sort your job when listing it
#SBATCH --cpus-per-task 8             # number of cpus for each task. One by default
#SBATCH --partition shared-cpu
#SBATCH --time 01:00:00
#SBATCH --mem-per-cpu=4000 # in MB
#SBATCH --output=slurm-output/%x-%A/%a.out
#SBATCH --array=0-365%10
#SBATCH --exclude cpu029,cpu030,cpu031,cpu032,cpu079,cpu123,cpu125,cpu137
PATH_TO_EXE=$HOME/dampe-stk-ml/inference.py
TMP_SSD_PATH_ON_NODE=/scratch/$USER/$SLURM_ARRAY_TASK_ID/
INPUT_DATA_PATH=$HOME/scratch/selectedData_Jennifer/
REGEX_DATA_FILE="DmlNtup_Data_2016[0-9][0-9]_photons[0-9][0-9]_MPhotonSelProj_TQ_SkimAllVar_FullSel\.root"
OUTPUT_DIR=$HOME/scratch/parsed_selected_validation/
LIST_OF_ALL_FILES="${TMP_SSD_PATH_ON_NODE}files_to_process.txt"
PATH_TO_CVMFS=/cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos8-gcc11-opt/setup.sh
COARSE_MODEL=$HOME/dampe-stk-ml/cluster/signal_only/inference/best_coarse/checkpoint-10-11.48.hdf5
COARSE_RESOLUTION=500 #um/pixel
FINE_MODEL=$HOME/dampe-stk-ml/cluster/signal_only/inference/best_fine/checkpoint-36-33.12.hdf5
FINE_RESOLUTION=50 #um/pixel
NUM_WORKERS=8
BATCH_SIZE=1024

#Load the modules needed to activate the cvmfs environment later on
module load GCCcore
module load Python
module load libreadline
#Print the task id and other info
echo "I am task_id " ${SLURM_ARRAY_TASK_ID} " on node " $(hostname)" in directory " $(pwd)

#Create the cache directory on the SSD of the cluster
srun -J "Cache folder" zsh -c "mkdir -p ${TMP_SSD_PATH_ON_NODE}"

#Generate a list of all files to process
#They are located in the INPUT_DATA_PATH and match the REGEX_DATA_FILE using bash
srun -J "Preparing" zsh -c "echo 'Starting with the creation of the file list' && \
        files=\$(ls ${INPUT_DATA_PATH} | grep -E '${REGEX_DATA_FILE}') && \
        if [[ -z \"\$files\" ]]; then \
                echo 'Error: No files matched the regular expression' >&2; \
                exit 1; \
        else \
                echo \"\$files\" > ${LIST_OF_ALL_FILES}; \
        fi"

#Launch the inference script on the task_id-th file in the list
srun -J "Parsing" zsh -c "echo 'Launching inference script' \
        && source \"${PATH_TO_CVMFS}\" \
        && python -u \"${PATH_TO_EXE}\" --root-file \"${INPUT_DATA_PATH}\"\"$(sed -n $((SLURM_ARRAY_TASK_ID + 1))p ${LIST_OF_ALL_FILES})\" \
         --coarse-model ${COARSE_MODEL} --coarse-resolution ${COARSE_RESOLUTION} --only-coarse-model  --num-workers ${NUM_WORKERS} --batch-size ${BATCH_SIZE} \
         --output-folder ${OUTPUT_DIR} "
```
## Expected Result
What did you expect to happen when running the steps above?
All the jobs to terminate successfully (as they did in the past).
## Actual Result
Job terminates with the error message:
```
9.out:9:srun: error: cpu088: task 0: Bus error (core dumped)
```


## Post 2 by @Genevieve.Savard (2024-01-12T18:09:00.591Z)

Hi,
Your script uses storage that’s located on scratch, so you can expect problems as the scratch storage partition is broken:
[2024] Current issues on HPC Cluster[[2024] Current issues on HPC Cluster](https://hpc-community.unige.ch//hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245) HPC ChangeLog[HPC ChangeLog](https://hpc-community.unige.ch/c/hpc-announce/hpc-changelog/9)
> Yggdrasil: Scratch storage issue
Dear users, 
there is an issue on scratch storage on Yggdrasil 

Our team is actively working on resolving this problem, and we apologize for any inconvenience caused. We will keep you updated on the progress.

Thank you for your understanding. 
Best regards, 
update 15.01.2024: 

We’ll stop the scratch storage completely until further notice.
We’ll send today the disks to the vendor and they’ll try to rebuild the RAID in their lab.
We’ll increase temporarily the…


## Post 3 by @parzival.nussbaum (2024-01-13T09:59:46.956Z)

That makes sense! Thanks!


## Post 4 by @Adrien.Albert (2024-01-17T15:06:07.368Z)

Hi @parzival.nussbaum[@parzival.nussbaum](https://hpc-community.unige.ch/u/parzival.nussbaum)
Genevieve.Savard:
> Your script uses storage that’s located on scratch, so you can expect problems as the scratch storage partition is broken:
I am not 100% sure the scratch issue the root cause. To be sure, could you try it on Baobab and let me know how it performs?
