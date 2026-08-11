# Matlab on Baobab with OpenOnDemand

- Source: https://hpc-community.unige.ch/t/3460

- Created: 2024-05-23T16:40:44.911Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Jade.Awada (2024-05-23T16:40:44.953Z)

Dear HPC team,
I’m contacting you because I’m having an issue using Matlab on Baobab with OpenonDemand.
All my files/scripts are located on my $HOME/scratch.
I chose the public iterative cpu partition, with 1 node and 6 cores, for 8 hours to run a Matlab script.
However, when I attempt to run it, it constantly crashes with the message: “session in a bad state”. The data file I am using is 1.7 GB, and the script performs several iterations with this data. The last session I tried was job ID 10209205 .
I have also tried using other partitions, but the webpage TurboVNC disconnects every time.
FYI, I also tried to run the same script reducing the data (47.8MB) and using the public iterative cpu partition, with 1 node and 4 cores, for 8 hours, and it worked perfectly well ( job ID 10208383 if it’s useful).
I guess the problem could come from the data size then? I am wondering if I am selecting the correct partition, number of cores, or the appropriate amount of memory per job.
Any guidance would be very appreciated.
Thank you in advance for your help!!
Best,
Jade


## Post 2 by @Adrien.Albert (2024-05-24T09:40:38.832Z)

Hi @Jade.Awada[@Jade.Awada](https://hpc-community.unige.ch/u/jade.awada)
You have exceeded your memory allocation (OUT_OF_MEMORY)
```
(baobab)-[root@login2 ~]$ sacct -X -o Jobid%15,jobname,account,user,nodelist,ReqCPUS,ReqMem,ntask,start,end,Elapsed,state%20  -j 10209205
          JobID    JobName    Account      User        NodeList  ReqCPUS     ReqMem   NTasks               Start                 End    Elapsed                State 
--------------- ---------- ---------- --------- --------------- -------- ---------- -------- ------------------- ------------------- ---------- -------------------- 
       10209205 sys/dashb+    sinanaj     awada          cpu007        4     12000M          2024-05-23T18:13:57 2024-05-23T18:21:51   00:07:54        OUT_OF_MEMORY 
```
You need to ask more memory (you can try with 16GB)
I invite you  to read this Best Pratice about resource allocations, to help you :
https://doc.eresearch.unige.ch/hpc/best_practices?s[]=resource#stop_wasting_resources[https://doc.eresearch.unige.ch/hpc/best_practices?s[]=resource#stop_wasting_resources](https://doc.eresearch.unige.ch/hpc/best_practices?s%5B%5D=resource#stop_wasting_resources)
---
PS: At some point, using the interactive session has its limits, and it is preferable to run jobs asynchronously with sbatch (slurm batch).
Example:
- Writte your Matlab script
```
(baobab)-[alberta@login2 Matlab_project]$ cat simple_script.m
% simple_script.m
n = 5;
result = sum((1:n).^2);
fileID = fopen('output.txt', 'w');
fprintf(fileID, 'The sum of the squares of the first %d integers is : %d\n', n, result);
fclose(fileID);
exit;
```
- Write your sbatch, specifying the necessary resources, loading the required modules, then run the command with slurm.
```
(baobab)-[alberta@login2 Matlab_project]$ cat simple_script.sbatch
#!/bin/bash
#SBATCH --job-name=simple_matlab_job
#SBATCH --output=simple_matlab_job.out
#SBATCH --error=simple_matlab_job.err
#SBATCH --time=00:05:00
#SBATCH --partition=debug-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

module load MATLAB
srun matlab -nodisplay -nosplash -nodesktop -r "simple_script; exit"
```
- Add your job in the queue:
```
(baobab)-[alberta@login2 Matlab_project]$ sbatch simple_script.sbatch
Submitted batch job 10235116
```
- wait the result drinking your favorite coffee (for me it’s roasted hazelnut latte)
Job pending → Preparing my Latte
```
(baobab)-[alberta@login2 Matlab_project]$ sac -j 10236031
          JobID    JobName    Account      User        NodeList   NTasks               Start                 End      State 
--------------- ---------- ---------- --------- --------------- -------- ------------------- ------------------- ---------- 
       10236031 simple_ma+      burgi   alberta   None assigned                      Unknown             Unknown    PENDING
```
Job Running  → Drinking my Latte
```
(baobab)-[alberta@login2 Matlab_project]$ sac -j 10236031
          JobID    JobName    Account      User        NodeList   NTasks               Start                 End      State 
--------------- ---------- ---------- --------- --------------- -------- ------------------- ------------------- ---------- 
       10236031 simple_ma+      burgi   alberta          cpu001          2024-05-24T11:32:48             Unknown    RUNNING
```
Job Completed → FInish my Latte and read my result
```
(baobab)-[alberta@login2 Matlab_project]$ sac -j 10236031
          JobID    JobName    Account      User        NodeList   NTasks               Start                 End      State 
--------------- ---------- ---------- --------- --------------- -------- ------------------- ------------------- ---------- 
       10236031 simple_ma+      burgi   alberta          cpu001          2024-05-24T11:32:48 2024-05-24T11:33:14  COMPLETED
```
```
(baobab)-[alberta@login2 Matlab_project]$ cat output.txt 
The sum of the squares of the first 5 integers is : 55
```
Best Regards


## Post 3 by @Jade.Awada (2024-05-31T11:46:21.252Z)

Hi @Adrien.Albert[@Adrien.Albert](https://hpc-community.unige.ch/u/Adrien.Albert),
Thank you for your reply and your help!
It’s all good then!!!
Best,
Jade
