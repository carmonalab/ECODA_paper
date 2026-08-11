# Baobab Application launch failed since yesterday afternoon

- Source: https://hpc-community.unige.ch/t/4345

- Created: 2026-07-10T07:02:41.843Z

- Tags: baobab, slurm

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Erica.Lastufka (2026-07-10T07:02:41.905Z)

Hello, since yesterday around 5 pm I’ve had all my slurm scripts fail to launch. They worked a few hours previous. I get error messages like:
```
srun: error: slurm_receive_msgs: \[gpu016.baobab:6818\] failed: Socket timed out on send/recv operation
srun: error: Task launch for StepId=10301666.2 failed on node gpu016: Socket timed out on send/recv operation
srun: error: Application launch failed: Socket timed out on send/recv operation
srun: Job step aborted
```
Additionally, I now have a few jobs stuck in “COMPLETING” status since then. Is something going on with the cluster?


## Post 2 by @Yann.Sagon (2026-07-10T15:28:25.992Z)

Hi,
How do you submit your jobs?
I tried right now without issue:
```
(base) (baobab)-[lastufka@login1 ~]$ srun hostname
cpu001.baobab
```


## Post 3 by @Erica.Lastufka (2026-07-10T15:43:36.368Z)

There are additional steps I have to do because my environments are often not used correctly by srun. I have this in my scripts:
unset PYTHONPATH
export PYTHONNOUSERSITE=1
export PATH=/home/users/l/lastufka/.conda/envs/dinov3/bin:$PATH
alias mysrun=‘srun /home/users/l/lastufka/.conda/envs/dinov3/bin/python3’
then I launch with mysrun
All my scripts with these commands worked fine until yesterday afternoon.
I still have 4 jobs with the state “COMPLETING” that show up when check squeue. Not sure if this is related.


## Post 4 by @Yann.Sagon (2026-07-13T07:34:06.863Z)

Hi,
Do you need an interactive session? If not, you should definitely use `sbatch` instead of `srun`. If there is an issue with login1, for example, or with the connection between login1 and your desktop, the job will be terminated.
You are using
Erica.Lastufka:
> unset PYTHONPATH
in your script. This variable is set because you have modules loaded. If you don’t need them, just use `module purge` to cleanly unload all your modules. For the record, you have defined the modules to be loaded when logging in in your `.bash_profile`. If you no longer need them, clean them up there too.


## Post 5 by @Erica.Lastufka (2026-07-13T08:56:15.434Z)

Thanks, I used module purge and commented out these lines:
unset PYTHONPATH
export PYTHONNOUSERSITE=1
Now the jobs launch and execute.
