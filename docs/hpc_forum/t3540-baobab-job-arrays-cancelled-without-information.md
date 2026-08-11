# Baobab, Job-arrays Cancelled without Information

- Source: https://hpc-community.unige.ch/t/3540

- Created: 2024-07-16T10:04:29.866Z

- Posts: 9

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Ludovic.Dumoulin (2024-07-16T10:04:29.905Z)

Hello,
Username: dumoulil
Cluster: Baobab
Subject: My jobs are getting cancelled
Job-array ID: 11389088, the 12 jobs that were running: 11389100, 11389099, 11389098, 11389097, 11389096, 11389095, 11389094, 11389093, 11389092, 11389091, 11389090, 11389089
I had to run 91 jobs, only 12 started and were cancelled, in the .out file I have the following message:
> slurmstepd: error: *** JOB 11389100 ON gpu031 CANCELLED AT 2024-07-15T20:01:00 ***
It is the second time that it happens with these jobs. Last time it was correponding tothe crash of:  /dpnc/beegfs/… But I am saving in /srv/beegfs/scratch/users/d/… is it related ?
Thank you for your help


## Post 2 by @Ludovic.Dumoulin (2024-07-24T12:43:47.490Z)

Hello,
The problem still persist, I have no idea to solve it…
> slurmstepd: error: *** JOB 11516164 ON gpu020 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516171 ON gpu030 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516172 ON gpu031 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516173 ON gpu031 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516177 ON gpu030 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516178 ON gpu030 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516379 ON gpu031 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516182 ON gpu020 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516184 ON gpu030 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516185 ON gpu031 CANCELLED AT 2024-07-24T12:01:00 ***
> slurmstepd: error: *** JOB 11516192 ON gpu020 CANCELLED AT 2024-07-24T12:01:00 ***
I find weird that the two arrays were cancelled at hh:01:00


## Post 3 by @Ludovic.Dumoulin (2024-07-24T13:37:37.150Z)

The issue is similar to this post, but no answer has been provided.
I ask him, the problem is still persisting…
Partial cancellation of jobs, others completed fine[Partial cancellation of jobs, others completed fine](https://hpc-community.unige.ch//hpc-community.unige.ch/t/partial-cancellation-of-jobs-others-completed-fine/3496) HPC issues[HPC issues](https://hpc-community.unige.ch/c/hpc-support/hpc-issues/14)
> Hi HPC team, 
Primary informations
Username: dedenon 
Cluster: Baobab 
Description
I ran a CPU job array last Friday on our private partition kruse-cpu, and I got cancellation of more than half of them, but the remaining ones are completed. 
 [jobarray][[jobarray]](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/b81e83470c2e3c1ad6bf78d23e5f3c9e31d07091.png) 
Here is the sbatch file 
 [sbatch][[sbatch]](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/6102f0df183a2fb8b04bd2d5ca9cc4e2ddee798d.png) 
I don’t think the issue is in my code, because I do parametric scanning for simulations with 20 independent realizations for 5 000 000 time steps, and some of them got completed (see slurm output below) 
job 9…
Is it related to rules of private kruse partitions ? But he is using CPU and I use GPU …


## Post 4 by @Adrien.Albert (2024-07-24T14:36:50.642Z)

Hi @Ludovic.Dumoulin[@Ludovic.Dumoulin](https://hpc-community.unige.ch/u/ludovic.dumoulin)
Could share your sbatch please ?


## Post 5 by @Ludovic.Dumoulin (2024-07-25T08:17:18.735Z)

Hello, Thanks, here is my sbatch
> #!/bin/env bash
> #SBATCH --array=1-91%40
> #SBATCH --partition=private-kruse-gpu
> #SBATCH --time=7-00:00:00
> #SBATCH --output=%J.out
> #SBATCH --mem=3000
> #SBATCH --gpus=ampere:1
> #SBATCH --constraint=DOUBLE_PRECISION_GPU
> module load Julia
> cd /srv/beegfs/scratch/users/d/dumoulil/Data/P-series/Dt/
> srun julia --optimize=3 /home/users/d/dumoulil/Code/FFT_2D_P_Dt/2D.jl
EDIT : In my previous, sbatch the there was also `shared-gpu` that I removed for this try (11531921). Is it possible that my jobs were cancelled because they were waiting in the shared partition for too long ?


## Post 6 by @Adrien.Albert (2024-07-30T10:32:51.526Z)

Hi @Ludovic.Dumoulin[@Ludovic.Dumoulin](https://hpc-community.unige.ch/u/ludovic.dumoulin)
```
(baobab)-[root@admin1 ~]$ sacct -X -o Jobid%15,jobname,account,user,nodelist,ReqCPUS,ReqMem,ntask,start,end,Elapsed,state%15 -j 11334268
          JobID    JobName    Account      User        NodeList  ReqCPUS     ReqMem   NTasks               Start                 End    Elapsed           State 
--------------- ---------- ---------- --------- --------------- -------- ---------- -------- ------------------- ------------------- ---------- --------------- 
     11334268_1      Dt.sh     krusek  dumoulil          gpu020        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
     11334268_2      Dt.sh     krusek  dumoulil          gpu020        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
     11334268_3      Dt.sh     krusek  dumoulil          gpu020        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
     11334268_4      Dt.sh     krusek  dumoulil          gpu030        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
     11334268_5      Dt.sh     krusek  dumoulil          gpu030        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
     11334268_6      Dt.sh     krusek  dumoulil          gpu030        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
     11334268_7      Dt.sh     krusek  dumoulil          gpu031        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
     11334268_8      Dt.sh     krusek  dumoulil          gpu031        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
     11334268_9      Dt.sh     krusek  dumoulil          gpu031        1      3000M          2024-07-12T15:35:19 2024-07-13T02:01:00   10:25:41  CANCELLED by 0 
    11334268_10      Dt.sh     krusek  dumoulil          gpu031        1      3000M          2024-07-12T15:38:46 2024-07-13T02:01:00   10:22:14  CANCELLED by 0 
    11334268_11      Dt.sh     krusek  dumoulil          gpu030        1      3000M          2024-07-12T15:38:49 2024-07-13T02:01:00   10:22:11  CANCELLED by 0 
    11334268_12      Dt.sh     krusek  dumoulil          gpu020        1      3000M          2024-07-12T15:38:51 2024-07-13T02:01:00   10:22:09  CANCELLED by 0 
11334268_[13-9+      Dt.sh     krusek  dumoulil   None assigned        1      3000M                         None 2024-07-13T02:01:00   00:00:00  CANCELLED by 0
```
```
$ grep 11334268 /var/log/slurm/*
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:18.183] _slurm_rpc_submit_batch_job: JobId=11334268 InitPrio=1560177 usec=603
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.004] sched: Allocate JobId=11334268_1(11334269) NodeList=gpu020 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.005] sched: Allocate JobId=11334268_2(11334270) NodeList=gpu020 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.006] sched: Allocate JobId=11334268_3(11334271) NodeList=gpu020 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.007] sched: Allocate JobId=11334268_4(11334272) NodeList=gpu030 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.007] sched: Allocate JobId=11334268_5(11334273) NodeList=gpu030 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.008] sched: Allocate JobId=11334268_6(11334274) NodeList=gpu030 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.009] sched: Allocate JobId=11334268_7(11334275) NodeList=gpu031 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.010] sched: Allocate JobId=11334268_8(11334276) NodeList=gpu031 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:35:19.011] sched: Allocate JobId=11334268_9(11334277) NodeList=gpu031 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:38:46.961] sched: Allocate JobId=11334268_10(11334287) NodeList=gpu031 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:38:49.448] sched: Allocate JobId=11334268_11(11334288) NodeList=gpu030 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-12T15:38:51.921] sched: Allocate JobId=11334268_12(11334289) NodeList=gpu020 #CPUs=1 Partition=private-kruse-gpu
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:01:00.083] _slurm_rpc_kill_job: REQUEST_KILL_JOB JobId=11334268 uid 0
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:02:32.526] cleanup_completing: JobId=11334268_12(11334289) completion process took 92 seconds
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:02:32.531] cleanup_completing: JobId=11334268_10(11334287) completion process took 92 seconds
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:02:32.533] cleanup_completing: JobId=11334268_4(11334272) completion process took 92 seconds
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:02:32.561] cleanup_completing: JobId=11334268_5(11334273) completion process took 92 seconds
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:02:32.572] cleanup_completing: JobId=11334268_2(11334270) completion process took 92 seconds
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:02:33.553] cleanup_completing: JobId=11334268_11(11334288) completion process took 93 seconds
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:02:33.560] cleanup_completing: JobId=11334268_9(11334277) completion process took 93 seconds
/var/log/slurm/slurmctld.log-20240715:[2024-07-13T02:02:33.724] cleanup_completing: JobId=11334268_3(11334271) completion process took 93 seconds
```
According to the logs and the precise TIMELIMIT, I suspect the jobs had a TimeLimit set to 10:25:00 and took too long to complete.
I thought the status for this kind of situation was “TIMELIMIT,” but perhaps “CANCELLED by 0” means “TIMELIMIT + too much time to complete”?
> completion process took 92 seconds
I do not find confirmation, but the completion process should take 90 seconds. After this time, Slurm kills the jobs to free up the node, which could explain the “Cancelled by 0” status.
The only thing I do not understand is your sbatch script, which sets the TIMELIMIT to 7 days. Are you sure you provided the exact sbatch script you executed for this job?
I have opened a ticket with Slurm to ensure this behavior is not associated with another event.
Best Regards


## Post 7 by @Ludovic.Dumoulin (2024-08-14T08:18:06.161Z)

Hello,
I am sorry my late reply, I was in vacation.
My exact bash was different because I forgot to remove the keyword `shared-gpu` from the partition list:
Ludovic.Dumoulin:
> was also `shared-gpu` that I removed for this try (11531921). Is it possible that my jobs were cancelled because they were waiting in the shared partition for too long ?
The old bash was:
```
#!/bin/env bash
#SBATCH --array=1-91%40
#SBATCH --partition=private-kruse-gpu,shared-gpu
#SBATCH --time=7-00:00:00
#SBATCH --output=%J.out
#SBATCH --mem=3000
#SBATCH --gpus=ampere:1
#SBATCH --constraint=DOUBLE_PRECISION_GPU
module load Julia
cd /srv/beegfs/scratch/users/d/dumoulil/Data/P-series/Dt/
srun julia --optimize=3 /home/users/d/dumoulil/Code/FFT_2D_P_Dt/2D.jl
```
Without `shared-gpu` the job array `11531921` is running fine.


## Post 8 by @Adrien.Albert (2024-08-14T12:15:27.071Z)

Hi @Ludovic.Dumoulin[@Ludovic.Dumoulin](https://hpc-community.unige.ch/u/ludovic.dumoulin)
We provide the log to Slurm support and they have no clue about this issue.
Has the problem happened again?


## Post 9 by @Ludovic.Dumoulin (2024-08-14T12:26:12.249Z)

Hello,
It hasn’t happened again. I’ll see for my next job array (after this one).
Thank you,
Best,
