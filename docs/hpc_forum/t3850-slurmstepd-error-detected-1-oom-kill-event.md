# Slurmstepd: error: Detected 1 oom_kill event

- Source: https://hpc-community.unige.ch/t/3850

- Created: 2025-03-05T14:27:07.219Z

- Tags: baobab

- Posts: 5

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Camille.Serquet (2025-03-05T14:27:07.261Z)

I am running a MATLAB code doing fMRI analysis. The input is a matrix of ~60GB that I need to load, and then multiple calculations are done.
The loading is ok.
The first few calculations are ok.
I clear all variables that are not useful any more.
All is stored on my scratch directory.
Working directory is also the scratch directory.
BUT, I still end up with this error :
> slurmstepd: error: Detected 1 oom_kill event in StepId=15279911.0. Some of the step tasks have been OOM Killed.
> srun: error: cpu245: task 0: Out Of Memory
Any idea what kills my job ? if not, is there a way to have more information about the kill ?
Here is my batch script (which worked for all previous MATLAB code I used)
```
#!/bin/bash

# Name of the job (optional I guess)
#SBATCH --job-name=kfind

# Send an email when the job is completed
#SBATCH --mail-user=camille.serquet@unige.ch
#SBATCH --mail-type=END

# Run the job at HH:MM today
#SBATCH --begin=now

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10

# Select partition, public-cpu for <4days, shared-cpu <12h, shared-bigmem, for big memory <12h
#SBATCH --partition=shared-bigmem

# Ask to use full memory of core when you have huge data
#SBATCH --mem=0

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=10:00:00

# Set the working directory - put your MATLAB script there
#SBATCH --chdir=/srv/beegfs/scratch/users/s/serquet/

# Request the MATLAB license 
#SBATCH --licenses=matlab@matlablm.unige.ch

# Load MATLAB on BAOBAB
module load MATLAB

# Set your MATLAB *.m file name
BASE_MFILE_NAME=mCAP_method2_findbestk_pca_Camille_all

# Run your MATLAB code and be happy 
srun matlab -nodesktop -nosplash -nodisplay -r ${BASE_MFILE_NAME}
```


## Post 2 by @Gael.Rossignol (2025-03-06T08:06:47.768Z)

Dear Camille,
I check your job with seff :
```
 (baobab)-[root@admin1 ~]$ seff 15279911
Job ID: 15279911
Cluster: baobab
User/Group: serquet/hpc_users
State: FAILED (exit code 1)
Nodes: 1
Cores per node: 10
CPU Utilized: 00:55:10
CPU Efficiency: 28.17% of 03:15:50 core-walltime
Job Wall-clock time: 00:19:35
Memory Utilized: 237.71 GB
Memory Efficiency: 95.08% of 250.00 GB (250.00 GB/node)
```
Then I see that you have an Out Of Memory issue, your job failed because you don’t have enough memory.
Multiple solutions are possible, I propose tu use nodes with more memory in the shared-cpu partition (cpu312,cpu313) as you job is less than 12h. Or migrate to bamboo where shared-bigmem partition has more memory on nodes.
Best regards,


## Post 3 by @Gael.Rossignol (2025-03-06T10:53:41.959Z)

Dear Camille,
As someone advise me, it’s preferable to fix the memory on the job and do not use :
```
#SBATCH --mem=0
```
Because job can be triggered on a node with poor memory ressources and fail.
Please update by using for exemple on your job :
```
#SBATCH --mem=400G
```
Best regards,


## Post 4 by @Camille.Serquet (2025-03-08T12:34:03.569Z)

Thanks @Gael.Rossignol[@Gael.Rossignol](https://hpc-community.unige.ch/u/Gael.Rossignol) !
The memory is sufficient when I set it to mem=400G !
BUT, previously, I set the wall time to 12:00:00 just because I didn’t know how much time it would actually take.
And now, I need it to be more than 12:00:00, with the same mem=400G.
Which cluster or/and partition will be best for the job ? (I tried a couple options on Baobab, but it keeps refusing my job :frowning:


## Post 5 by @Gael.Rossignol (2025-03-10T08:37:16.128Z)

Dear Camille,
If you need more than 12h of processing and high amount of memory, the only way is to migrate on bamboo and use public-bigmen :slight_smile:
```
PartitionName=public-gpu
   AllowGroups=ALL AllowAccounts=ALL AllowQos=ALL
   AllocNodes=ALL Default=NO QoS=N/A
   DefaultTime=00:02:00 DisableRootJobs=NO ExclusiveUser=NO ExclusiveTopo=NO GraceTime=0 Hidden=NO
   MaxNodes=UNLIMITED MaxTime=2-00:00:00 MinNodes=0 LLN=NO MaxCPUsPerNode=UNLIMITED MaxCPUsPerSocket=UNLIMITED
   Nodes=gpu[001-003]
   PriorityJobFactor=1 PriorityTier=1 RootOnly=NO ReqResv=NO OverSubscribe=NO
   OverTimeLimit=NONE PreemptMode=OFF
   State=UP TotalCPUs=272 TotalNodes=3 SelectTypeParameters=NONE
   JobDefaults=(null)
   DefMemPerNode=UNLIMITED MaxMemPerNode=UNLIMITED
   TRES=cpu=266,mem=1500G,node=3,billing=773,gres/gpu=20,gres/gpu:nvidia_a100_80gb_pcie=4,gres/gpu:nvidia_geforce_rtx_3090=16
   TRESBillingWeights=CPU=1.0,Mem=0.25G,GRES/gpu=1,GRES/gpu:nvidia_a100-pcie-40gb=5,GRES/gpu:nvidia_a100_80gb_pcie=8,GRES/gpu:nvidia_geforce_rtx_2080_ti=2,GRES/gpu:nvidia_geforce_rtx_3080=3,GRES/gpu:nvidia_geforce_rtx_3090=5,GRES/gpu:nvidia_geforce_rtx_4090=8,GRES/gpu:nvidia_rtx_a5000=5,GRES/gpu:nvidia_rtx_a5500=5,GRES/gpu:nvidia_rtx_a6000=8,GRES/gpu:nvidia_titan_x=1,GRES/gpu:tesla_p100-pcie-12gb=1
   ResumeTimeout=GLOBAL SuspendTimeout=GLOBAL SuspendTime=GLOBAL PowerDownOnIdle=NO
```
This allow 2 days of processing on nodes with 1T of memory.
Best regards,
