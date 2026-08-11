# How to improve matlab-based codes in HPC UNIGE

- Source: https://hpc-community.unige.ch/t/3910

- Created: 2025-04-04T11:03:22.299Z

- Posts: 2

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Julian.GaviriaLopez (2025-04-04T11:03:22.342Z)

I ran some MATLAB clustering analyses in Yggdrasil. Parallel computing in MATLAB is highly recommended to improve the performance. I first compared my code between parallel computing (threads vs. process ( see ref here)[here)](https://ch.mathworks.com/help/parallel-computing/quick-start-parallel-computing-in-matlab.html). The results yielded faster performance in parallel computing with threads. Apparently, because of some overhead of copying data to the workers (i.e., the cores).
Then I compared the same run in the following scenario:  Parallel computing (threads) vs. no parallel computing. Again, the results yielded better performance of parallel computing threads in terms of time. A last comparison with two different parallel comp threads (5 cores vs. 18 cores) also indicated that more cores increase the efficiency of the code, specially in the elapsed time.  The results of the mentined comparisons are described below.
Questions:
- Can I request 16 cores for the big-mem partition? I’ve tried to run my code with the final sample size (126 Gb) with 2 cores. Unfortunately, 4 days was not enough.
- is there any alternative to monitor CPU and RAM usage different than the “seff” command?  “sstat” did not work, and “sacct” yield outputs difficult to interpret.
Thanks in advance for your comments.
```
Cores used in MATLAB: 1 (default)
Elapsed time is 60.777525 seconds.
Job ID: 39752917
Cluster: yggdrasil
User/Group: gavirial/unige
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:06:08
CPU Efficiency: 4.42% of 02:18:40 core-walltime
Job Wall-clock time: 00:06:56
Memory Utilized: 7.94 GB
Memory Efficiency: 13.23% of 60.00 GB (60.00 GB/node)
```
---
```
Cores used in MATLAB: 18
Job ID: 39750731
Elapsed time is 17.334834 seconds.
Cluster: yggdrasil
User/Group: gavirial/unige
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 00:05:30
CPU Efficiency: 6.60% of 01:23:20 core-walltime
Job Wall-clock time: 00:04:10
Memory Utilized: 13.71 GB
Memory Efficiency: 22.84% of 60.00 GB (60.00 GB/node)
```
---
```
Cores used in MATLAB: 5
Elapsed time is 23.345938 seconds.
Job ID: 39753389
Cluster: yggdrasil
User/Group: gavirial/unige
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 30
CPU Utilized: 00:08:55
CPU Efficiency: 6.07% of 02:27:00 core-walltime
Job Wall-clock time: 00:04:54
Memory Utilized: 7.34 GB
Memory Efficiency: 12.23% of 60.00 GB (60.00 GB/node)
```


## Post 2 by @Yann.Sagon (2025-04-08T16:22:04.350Z)

Dear @Julian.GaviriaLopez[@Julian.GaviriaLopez](https://hpc-community.unige.ch/u/julian.gavirialopez),
Julian.GaviriaLopez:
> A last comparison with two different parallel comp threads (5 cores vs. 18 cores) also indicated that more cores increase the efficiency of the code, specially in the elapsed time.
If I understand correctly: you launched the same job three time with the following parameters:
- 1 core 60sec
- 5 cores 23sec
- 18 cores 17sec
Thanks chatgpt for the help with this graph:
e76676d1-e5f5-467b-9c7a-b308417adf90
e76676d1-e5f5-467b-9c7a-b308417adf902779×980 130 KB
[e76676d1-e5f5-467b-9c7a-b308417adf902779×980 130 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/ed32548ca60d3ac65035ae8da463f984b4bafcf3.png)
You can see that the efficiency is dropping by 50% when using 5 cores: it means the code you are using isn’t very well parallelized and you shouldn’t use more than 5 cores. To be exhaustive you should re run your test with 2 and 4 cores to see which number is more efficient.
In the seff’s output, it is indicated you allocated 20 and 30 cores: you should only allocate the number you’ll use with matlab (i.e. 5 for example), and if you want to do a fair comparison between runs, you should ideally use the same compute node generation.
Julian.GaviriaLopez:
> Can I request 16 cores for the big-mem partition? I’ve tried to run my code with the final sample size (126 Gb) with 2 cores
As all our newer compute nodes have 128 cores and 512GB of RAM, you can simply request resource in the partition public-cpu with 4 cores and the required RAM (on Bamboo).
Julian.GaviriaLopez:
> is there any alternative to monitor CPU and RAM usage different than the “seff” command? “sstat” did not work, and “sacct” yield outputs difficult to interpret.
`sstat` is working, but you need to specify the job step or option `--all` to display the information for every job step of your job.
Example:
```
(yggdrasil)-[root@admin1 ~]$ sstat -j 39919765 --all | less -#2 -N -S
      1 JobID         MaxVMSize  MaxVMSizeNode  MaxVMSizeTask  AveVMSize     MaxRSS MaxRSSNode MaxRSSTask     AveRSS MaxPages MaxPagesNode   MaxPagesTask   AvePages     MinCPU MinCPUNode MinCPUTask     AveCPU   NTasks AveCPUFreq ReqCPUF>
      2 ------------ ---------- -------------- -------------- ---------- ---------- ---------- ---------- ---------- -------- ------------ -------------- ---------- ---------- ---------- ---------- ---------- -------- ---------- ------->
      3 39919765.ex+    489808K         cpu095              0      5584K      1728K     cpu095          0      1728K        0       cpu095              0          0   00:00:00     cpu095          0   00:00:00        1      2.80M        >
      4 39919765.ba+   1298188K         cpu095              0   1298188K    599588K     cpu095          0    599588K     2511       cpu095              0       2511   00:39:25     cpu095          0   00:39:25        1       984K        >
```
What is the issue with `seff`?
