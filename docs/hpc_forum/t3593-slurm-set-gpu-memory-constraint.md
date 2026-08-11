# SLURM set GPU memory constraint

- Source: https://hpc-community.unige.ch/t/3593

- Created: 2024-08-14T12:05:45.168Z

- Tags: bamboo

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @daniel.forerosanchez (2024-08-14T12:05:45.229Z)

Hi, I’ve been trying to test the large memory GPU on bamboo using
`srun -n1 -c32 -p public-gpu --mem-per-gpu=80G --gpus=1 --gres=gpu:1 -t 1:00:00 python ...`
and while it works, I keep getting a resource exhausted message from JAX when it tries to allocate ~15Gb, how can I make sure I am using the right node and perhaps check my total memory usage? I’ve checked and it usually allocates to gpu001 node which does not have the large memory GPUs. Is there any way I can just test if my code will run there? (i.e. if memory is enough for my problem)?
Thanks in advance


## Post 2 by @Adrien.Albert (2024-08-14T14:56:30.345Z)

Hi Daniel,
When you start a job with Slurm, it gets saved in the `slurmdbd` . You can check the status or details of your job using commands like `sacct` or `squeue`:
- `squeue`: Shows information about jobs in the Slurm scheduling queue.
- `sacct`: Displays detailed accounting data for jobs in the Slurm job accounting log or database.
For instance, to get detailed information about your job (59075), you can use:
```
sacct -X -o Jobid%15,jobname,account,user,nodelist,ReqCPUS,ReqMem,ntask,start,end,Elapsed,state -j 59075
         JobID    JobName    Account      User        NodeList  ReqCPUS     ReqMem   NTasks               Start                 End    Elapsed      State 
--------------- ---------- ---------- --------- --------------- -------- ---------- -------- ------------------- ------------------- ---------- ---------- 
         59075     python       pepe  dforeros          gpu003       10     30000M          2024-08-14T14:07:02 2024-08-14T14:10:13   00:03:11  COMPLETED
```
You can also use the seff command, which gives a summary of the job’s efficiency:
```
(bamboo)-[root@admin1 ~]$ seff 59075
Job ID: 59075
Cluster: bamboo
User/Group: dforeros/hpc_users
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 10
CPU Utilized: 00:01:08
CPU Efficiency: 3.56% of 00:31:50 core-walltime
Job Wall-clock time: 00:03:11
Memory Utilized: 25.31 GB
Memory Efficiency: 86.40% of 29.30 GB
```
Just a heads-up, though: seff reports job stats after the job completes, so it gives a summary rather than real-time data.
What does mean large memory GPU :question:*
Here the gpu specification on Bamboo:
Node
Type
Model
GPU Memory (VRam)
gpu count
gpu[001-002]
Ampere
NVIDIA GeForce RTX 3090
24 GB
8
gpu003
Ampere
NVIDIA A100 80GB PCIe
80GB
4
It’s also important to differentiate between GPU memory (VRAM) and system memory (RAM).
If you specify `--gres=gpu:1,VramPerGpu:80G` you will exclude  RTX 3090 node.
As you can see there are several types of GPU, including one with “large memory” (:question:*) A100;
But before opting for it, it’s important to ask yourself whether you really need it. The A100 is about 20 times more expensive than the RTX 3090, and there are fewer of them in the cluster. So, you might end up blocking a valuable resource without gaining much in performance, depending on your workload.
For example, in job 59075, seff shows that you’re using almost 25GB of system memory. You might find that running the job on gpu[001-002] (with the RTX 3090) won’t significantly impact performance. It’s something you could test:
```
srun --ntasks=1 --cpus-per-task=1 --nodes=1 --nodelist=gpu[001-002] --partition=public-gpu --mem=25G --gres=gpu:1,VramPerGpu:24G Python ...
```
BUT it’s clear that we’d rather have an operational server than a standby server for certain jobs requiring only A100s. As I don’t have a precise idea of the number of jobs requiring this prerogative, I can’t say. So I’d say it depends on your patience and the urgency of your project.
(I have to admit that my answer is nuanced and raises other questions.)
