# Nested sbatch files submission through R

- Source: https://hpc-community.unige.ch/t/3499

- Created: 2024-06-18T13:35:32.599Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Davide.Pietrobon (2024-06-18T13:35:32.652Z)

Dear all,
I am currently working on a project where I need to execute multiple nested sbatch jobs from an R script on the Yggdrasil cluster. The R script (`test_ja_SA.R`) contains the following line of code:
`system("sbatch job_test_ja_wrapper.sh", intern = TRUE)`
Running the `sbatch job_test_ja_wrapper.sh` command directly from Putty successfully completes the job. The wrapper script, as suggested by its name, orchestrates the execution of several other sbatch files. All nested jobs (called `job_test_ja_row_count.sh` and `job_test_ja.sh`) run as expected. Please, find below the three sbatch files:
`sbatch job_test_ja_wrapper.sh`:
```
#!/bin/sh
#SBATCH --job-name=wrapper
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=1000
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=shared-cpu
# #SBATCH --mail-user=davide.pietrobon@unige.ch
# #SBATCH --mail-type=ALL

# Submit the row count job and wait for it to complete
JOBID_COUNT=$(sbatch job_test_ja_row_count.sh | cut -d ' ' -f 4)
echo "Row count job submitted with Job ID: $JOBID_COUNT"

# Wait for the row count job to finish
while squeue | grep -q $JOBID_COUNT; do
  sleep 10
done

# Read the number of rows
NUM_ROWS=$(cat row_count.txt)
echo "Number of rows: $NUM_ROWS"

# Submit the second job with the correct job array size and capture its Job ID
JOBID_ARRAY=$(sbatch --array=1-$NUM_ROWS job_test_ja.sh | cut -d ' ' -f 4)
echo "Array job submitted with Job ID: $JOBID_ARRAY"

echo $JOBID_ARRAY > job_array_id.txt

# Add a dependency to wait for the job array to complete successfully
sbatch --dependency=afterok:$JOBID_ARRAY --wrap="echo 'Job array $JOBID_ARRAY completed successfully'"`
```
`job_test_ja_row_count.sh`:
```
#!/bin/sh
#SBATCH --job-name=row_count
#SBATCH --time=00:05:00  # Set a shorter time as this job should be quick
#SBATCH --mem-per-cpu=1000  # Less memory might be required for this task
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=shared-cpu  # Specify the appropriate partition

module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1

# Run the R script that writes the number of rows to a file

INFILE="test_ja_row_count.R"

srun R CMD BATCH --no-save --no-restore "$INFILE"
```
`job_test_ja.sh`:
```
#!/bin/sh

#SBATCH --job-name=test_ja
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=5000        # Memory per CPU in MB
#SBATCH --cpus-per-task=1   # Each task gets 1 core, adjust this as needed
#SBATCH --ntasks=1        # Number of tasks you want to use
#SBATCH --partition=shared-cpu
#SBATCH --mail-user=davide.pietrobon@unige.ch
#SBATCH --mail-type=ALL

module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1

# Setup environment variable to pass to R script
export SLURM_ARRAY_TASK_ID

INFILE="test_ja.R"
# OUTFILE="test_ja_${SLURM_ARRAY_TASK_ID}.out"  # Output file includes task ID

# srun R CMD BATCH --no-save --no-restore "$INFILE" "$OUTFILE"
 srun R CMD BATCH --no-save --no-restore "$INFILE"
```
However, executing the same sbatch file from within R using `system()` results in failures specifically related to the job array managed by `job_test_ja_wrapper.sh`. The output suggests issues related to CPU binding (more on this below), which do not occur when the script is run from Putty. The sbatch file I use to run the R code is called `job_test_ja_SA.sh`:
```
#!/bin/sh
#SBATCH --job-name=ja_SA
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=shared-cpu
#SBATCH --mail-user=davide.pietrobon@unige.ch
#SBATCH --mail-type=ALL

module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1

# Run the R script that writes the number of rows to a file

INFILE="test_ja_SA.R"
OUTFILE="test_ja_SA.out"

srun R CMD BATCH --no-save --no-restore "$INFILE" "$OUTFILE"
```
Here is some of the output files:
- `slurm-34005396.out` refers to the output from running the `sbatch job_test_ja_wrapper.sh` command directly from Putty.
- `slurm-34005297.out` refers to the output from running the `sbatch job_test_ja_SA.sh`, which is the sbatch file that runs the R script `test_ja_SA.R` (i.e., the main R code that contains the line `system("sbatch job_test_ja_wrapper.sh", intern = TRUE)`).
- `slurm-34005308_15.out` is the output from an instance (task number 15—the output of each task is the same) of the job array as run by the `sbatch job_test_ja_wrapper.sh`.
`slurm-34005396.out`:
```
Job array 34005395 completed successfully
```
`slurm-34005297.out`:
```
srun: error: cpu119: task 0: Exited with exit code 1
```
`slurm-34005308_15.out`:
```
srun: error: CPU binding outside of job step allocation, allocated CPUs are: 0x00000000000000000010000000000000.
srun: error: Task launch for StepId=34005325.0 failed on node cpu124: Unable to satisfy cpu bind request
srun: error: Application launch failed: Unable to satisfy cpu bind request
srun: Job step aborted
```
Since this is my first experience with running nested sbatch jobs, I suspect I might be missing a crucial step or configuration. Any guidance or insights you could provide would be greatly appreciated.
Thank you so much for your kind support,
Davide


## Post 2 by @Davide.Pietrobon (2024-06-18T13:54:36.727Z)

Could the issue stem from submitting jobs to the queue from a compute resource? This possibility is suggested in discussions like the one found here: https://groups.google.com/g/slurm-users/c/mp_JRutKmCc?pli=1[https://groups.google.com/g/slurm-users/c/mp_JRutKmCc?pli=1](https://groups.google.com/g/slurm-users/c/mp_JRutKmCc?pli=1). Notably, submitting `sbatch job_test_ja_wrapper.sh` from a login node appears to work without errors, whereas attempts to submit it via an R script lead to problems.


## Post 3 by @Adrien.Albert (2024-06-18T14:50:33.815Z)

Hi Davide,
Following the google group post, did you tried the solution suggested here:
support.schedmd.com[support.schedmd.com](https://support.schedmd.com/show_bug.cgi?id=14298)
### 14298 – CPU binding outside of job step allocation[14298 – CPU binding outside of job step allocation](https://support.schedmd.com/show_bug.cgi?id=14298)
SchedMD - Slurm development and
                                 support. Providing support for some
                                 of the largest clusters in the
                                 world.
> Chris: also seen this recently under 22.05.  I think the issue is SLURM_CPU_BIND being inherited when sbatch is invoked and there therefore sometimes being a mismatch between the value of SLURM_CPU_BIND in the batch job and the taskset of the batch job: if you ‘unset SLURM_CPU_BIND’ before running sbatch then the issue doesn’t seem to occur.
> It seems like this is a change in behaviour in 22.05, but I’m not sure what’s caused it. Possibly a side effect of one of the following changes:
> – Fail srun when using invalid --cpu-bind options (e.g. --cpu-bind=map_cpu:99
> when only 10 cpus are allocated).
> – srun --overlap now allows the step to share all resources (CPUs, memory, and
> GRES), where previously --overlap only allowed the step to share CPUs with
