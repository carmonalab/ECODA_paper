# Job failed without .out file

- Source: https://hpc-community.unige.ch/t/3655

- Created: 2024-09-30T07:59:14.687Z

- Tags: baobab

- Posts: 3

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Xin.Wen (2024-09-30T07:59:14.729Z)

## Primary informations
Username: $wenx
Cluster: $baobab
## Description
Hello,
Hope you are doing well!
I have an issue as follows when running the jobs on HPC.
jobid: 12792706
Detail:
I was trying running the following .sh file
(Path:/home/users/w/wenx/stones/stones_forward_looking)
```
#!/bin/bash*
#SBATCH --partition=shared-cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=1000 # in MB
#SBATCH -o run_stones_sa-%A_%a.out
# Activate environment
ml GCCcore/11.2.0
ml Gurobi/9.5.0
ml Python/3.9.6
. ~/baobab_python_env/bin/activate
# Send costoptimal run to cluster
echo "Running Stones SA on node " $(hostname)
~/baobab_python_env/bin/python 02_sensitivity_analysis.py
```
However, I cannot find the .out file. When I check the job by running $seff 12792706:
```
(baobab)-[wenx@login1 stones_forward_looking]$ seff 12792706
perl: warning: Setting locale failed.
perl: warning: Please check that your locale settings:
LANGUAGE = (unset),
LC_ALL = (unset),
LC_CTYPE = "UTF-8",
LANG = "en_US.UTF-8"
are supported and installed on your system.
perl: warning: Falling back to a fallback locale ("en_US.UTF-8").
Job ID: 12792706
Cluster: baobab
User/Group: wenx/hpc_users
State: FAILED (exit code 15)
Nodes: 1
Cores per node: 4
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:00:00 core-walltime
Job Wall-clock time: 00:00:00
Memory Utilized: 1.03 MB
Memory Efficiency: 0.10% of 1000.00 MB
```
Also since I cannot find the .out file, I don’t know what went wrong.
Could you please help me with the issue and let me know why the .out file is not generated?
Please let me know if you need more information.
Thank you in advance for your help!
Best regards,
Xin


## Post 2 by @Yann.Sagon (2024-09-30T12:22:07.812Z)

Dear @Xin.Wen[@Xin.Wen](https://hpc-community.unige.ch/u/xin.wen)
the issue you are facing was due to an error with the storage. [2024] Current issues on HPC Cluster - #25 by Yann.Sagon[[2024] Current issues on HPC Cluster - #25 by Yann.Sagon](https://hpc-community.unige.ch/t/2024-current-issues-on-hpc-cluster/3245/25)
Sorry for the inconvenience, just resubmit your job. By the way: you are requesting only 1GB for your job, is it enough? By default, it is 3GB per core, thus you are requesting less than the default (4 x 3 => 12GB).
Best


## Post 3 by @Xin.Wen (2024-10-01T07:39:01.093Z)

Dear Yann,
Thank you so much for your prompt reply! It works now :slight_smile:  And Thank you for your advice!
Best,
Xin
