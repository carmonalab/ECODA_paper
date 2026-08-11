# Use Beast checkpointing with slurm job array

- Source: https://hpc-community.unige.ch/t/3771

- Created: 2024-12-19T10:40:37.653Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Cody.Cardenas (2024-12-19T10:40:37.653Z)

I found that pretty quickly on the HPC docs, but thanks for sharing it here.
I have started spinning my wheels on submitting a job that restarts until its completed. But I am getting a little stuck. I hope you, or someone else can help me with it. Let me know if I should start a new thread for it.
The checkpoints available for beast[checkpoints available for beast](https://beast.community/online_inference#loading-a-previously-generated-state-or-checkpoint-file-in-beast) seem pretty straight forward:
`beast  -load_state filename  -save_state newfilename [other commands] input.xml`
But I need to tell my script to look at the old `load_state` to make a new `save_state`; because if the log file exists beast wont start (`input.xml` turns into `input.log`).
my original slurm script is a pretty simple array that submits 16 jobs for 8 independent beast runs for two parallel analyses:
```
#!/bin/bash
#SBATCH --time 96:00:00 # in hours
#SBATCH --partition public-cpu
#SBATCH --cpus-per-task 16 # number of cores to use
#SBATCH --mem 80000 # in MB
#SBATCH --array 1-16 # run N jobs
#SBATCH --job-name BEAST_rerun

# to run: sbatch submit_jobs.sh
# be sure to adjust the length of your array based on your list

# get yer modules loaded
module load GCC/10.2.0 beagle-lib/3.1.2 Beast/1.10.5pre

# set variables
WORKING_DIR="/home/users/c/cardenac/phylo/2023_adephaga/BEAST_rerun"
ALIGN_LIST="jobs.list"

# your list should look like this:

# run01/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# run01/core_SORTADATE_beast_lognormal_sans-quasicalathus.xml
# run02/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# run02/core_SORTADATE_beast_lognormal_sans-quasicalathus.xml
# run03/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# ...

# get info from list & use awk to print that line based on the array where the line is based on the ${SLURM_ARRAY_TASK_ID}
DATASET_DIR=$(cat ${ALIGN_LIST} | awk -v var=${SLURM_ARRAY_TASK_ID} 'NR==var {print $0}' | cut -d "/" -f 1)
XML_FILE=$(cat ${ALIGN_LIST} | awk -v var=${SLURM_ARRAY_TASK_ID} 'NR==var {print $0}' | cut -d "/" -f 2)

# change directories to where data lives and run script
cd ${WORKING_DIR}/${DATASET_DIR}
beast -threads 16 ${XML_FILE}
```
I thought maybe using some kind of logical statement for writing a name to the load_state and file_state inputs might work and then wrapping my original command in srun… but I’m starting to get a little lost
```
#SBATCH --time 04:00:00 # in hours
#SBATCH --partition public-cpu
#SBATCH --cpus-per-task 16 # number of cores to use
#SBATCH --mem 80000 # in MB
#SBATCH --job-name BEAST_rerun
#SBATCH --array 1-16 # run N jobs

# to run: sbatch submit_jobs.sh
# be sure to adjust the length of your array based on the number of lines in your list
module load GCC/10.2.0 beagle-lib/3.1.2 Beast/1.10.5pre
# set variables
WORKING_DIR="/home/users/c/cardenac/phylo/2023_adephaga/BEAST_rerun"
ALIGN_LIST="jobs.list"

# your jobs.list should look like this:
# run01/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# run01/core_SORTADATE_beast_lognormal_sans-quasicalathus.xml
# run02/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# run02/core_SORTADATE_beast_lognormal_sans-quasicalathus.xml
# run03/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# ...

# get alignment path from list
DATASET_DIR=$(cat ${ALIGN_LIST} | awk -v var=${SLURM_ARRAY_TASK_ID} 'NR==var {print $0}' | cut -d "/" -f 1)
# use awk to print that line based on the array, where the line is based on the ${SLURM_ARRAY_TASK_ID}
XML_FILE=$(cat ${ALIGN_LIST} | awk -v var=${SLURM_ARRAY_TASK_ID} 'NR==var {print $0}' | cut -d "/" -f 2)
# set your prefix for checkpointing
PREFIX=$(cat ${ALIGN_LIST} | awk -v var=${SLURM_ARRAY_TASK_ID} 'NR==var {print $0}' | cut -d "/" -f 2 | cut -d "." -f 2 )

# need to append the job step to the run otherwise beast wont restart due to "the log file already exists"
if [ ! -f finished_${PREFIX} ] ; then
    # create a file to count the number of times beast has run
    # first check that the count file has been made
    if [ ! -f ${XML_FILE}_jobcount ] ; then
        touch ${XML_FILE}_jobcount.list
    fi

    # add your first XML file to your jobcount.list
	ls -1 ${XML_FILE} >> ${XML_FILE}_jobcount.list

    # counts the number of lines in your jobcount.list and makes a a count prefix
    COUNT_PREFIX=$(awk '{LIST_COUNT=FNR+1} {print LIST_COUNT}' ${XML_FILE}_jobcount.list)

    # okay neat we've got counts for a prefix, figure out how to call the original file
    COUNT_RUNS= #hmmm

    # maybe use the wrap function with srun?
    srun --dependency=afterany:$SLURM_JOBID \
#     --wrap ' # silencing this for my text editor!!
    cd ${WORKING_DIR}/${DATASET_DIR}
    beast -threads 16 -prefix ${XML_FILE}_run${COUNT_RUNS} -load_state filename[????] -save_state ${XML_FILE}_run${COUNT_RUNS}_checkpoint
#    '
fi
```
Am I over thinking it?
admin edit: created as new post


## Post 2 by @Cody.Cardenas (2024-12-19T15:30:59.019Z)

### Summary of previous post
I cant figure out how to submit a dependency[dependency](https://slurm.schedmd.com/sbatch.html#OPT_afternotok:job_id%5B:jobid...%5D) for an array job with the software beast[beast](https://beast.community/online_inference#loading-a-previously-generated-state-or-checkpoint-file-in-beast).
e.g. : `#SLURM --dependency=afternotok:$SLURM_JOBID` with
Because beast requires you to rename the checkpoint file with `-load_state filename -save_state filename_checkpoint` files. Because I need to start this 20+ times with a short run time (say 6-12 hours) I would need to inform the slurm `--dependency` command to look at a different `load_state` and `save_state` file for the next run: e.g,
- first run `... -load_state checkpoint_run1  -save_state checkpoint_run2...`
- second run `... -load_state checkpoint_run2  -save_state checkpoint_run3...`
- nth run `... -load_state checkpoint_run#N  -save_state checkpoint_run#N...`
But its not immediately clear how to do that; I’ve provided my attempt (see previous post) but I can’t wrap my head around how to tell slurm look at a different file, but keep restarting this job until it gives an OK exit status. I’m hoping there is a clearer way to do this rather than submitting 20-30 jobs for 16 independent runs.
---
### For now:
Since I kind of need to get this started before the end of the year, and there is compute space free on Bamboo, I’m gonna be very hacky about it and just monitor my jobs… maybe it’ll be helpful for someone if they have a better solution.
```
#!/bin/bash
#SBATCH --time 36:00:00 # in hours
#SBATCH --partition public-cpu
#SBATCH --cpus-per-task 96 # number of cores to use
#SBATCH --mem 80000 # in MB
#SBATCH --job-name BEAST_rerun_1
#SBATCH --array 1-16 # run N jobs

# cant figure out how to monitor a job with --dependency=afternotok:$SLURM_JOBID
# first run: sbatch submit_jobs.sh 01
# next runs: sbatch submit_jobs.sh 02 (see comment after first if-else statement)
# be sure to adjust the length of your array based on the number of lines in your list
module load GCC/10.2.0 beagle-lib/3.1.2 Beast/1.10.5pre

# set variables
WORKING_DIR="/home/users/c/cardenac/phylo/2023_adephaga/BEAST_rerun"
JOB_LIST="jobs.list"
# your jobs.list should look like this:
# run01/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# run01/core_SORTADATE_beast_lognormal_sans-quasicalathus.xml
# run02/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# run02/core_SORTADATE_beast_lognormal_sans-quasicalathus.xml
# run03/core_SORTADATE_beast_exponential_sans-quasicalathus.xml
# ...

# get alignment path from list
DATASET_DIR=$(cat ${JOB_LIST} | awk -v var=${SLURM_ARRAY_TASK_ID} 'NR==var {print $0}' | cut -d "/" -f 1)
# use awk to print that line based on the array, where the line is based on the ${SLURM_ARRAY_TASK_ID}
XML_FILE=$(cat ${JOB_LIST} | awk -v var=${SLURM_ARRAY_TASK_ID} 'NR==var {print $0}' | cut -d "/" -f 2)
# set your prefix for checkpointing
PREFIX=$(cat ${JOB_LIST} | awk -v var=${SLURM_ARRAY_TASK_ID} 'NR==var {print $0}' | cut -d "/" -f 2 | cut -d "." -f 2)
#take user input for sbatch script
LOAD_STATE=${1}
SAVE_STATE=$(echo ${LOAD_STATE} | awk -v PREVRUN=${LOAD_STATE} '{COUNT=PREVRUN+1} {print COUNT}')

# change into working directory
cd ${WORKING_DIR}/${DATASET_DIR}

beast -threads 96 \
	-prefix ${PREFIX}_run1 \
	-save_state ${PREFIX}_run1_checkpoint \
	-save_every 500000 \
	${XML_FILE}

# once your job has finished comment out your if-else step and run with this comented out section
#beast -threads 96 -prefix ${XML_FILE}_run${CUR_RUN} -load_state ${XML_FILE}_run${PREV_RUN}_checkpoint -save_state ${XML_FILE}_run${CUR_RUN}_checkpoint ${XML_FILE}
```
