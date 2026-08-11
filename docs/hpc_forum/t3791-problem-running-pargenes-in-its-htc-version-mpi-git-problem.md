# Problem running pargenes in its htc version (MPI/git problem)

- Source: https://hpc-community.unige.ch/t/3791

- Created: 2025-01-21T16:00:14.588Z

- Posts: 3

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Camille.Christe (2025-01-21T16:00:14.635Z)

Hi,
I try to use a ParGenes, a program that allow its usage in hpc.
It is installed via conda.
When running it on a single thread, it works fine but very slowly (pargenes.py). So I tried the hpc version (pargenes-hpc.py).
It is specified that it needs MPI to be installed. What I did with `ml GCC/13.2.0  OpenMPI/4.1.6`
It get the following error message, related to git.
This is the error message for job # 37669149
```
[Error] [0:00:00] mpi-scheduler execution failed with error code 1
[Error] [0:00:00] Will now exit...
[Error] <class 'RuntimeError'> mpi-scheduler  execution failed with error code 1
Writing report file in /home/christec/Clematis/Run3/Analyses/Tree/raxml_all/report.txt
When reporting the issue, please always send us this file.
fatal: not a git repository (or any parent up to mount point /)
Stopping at filesystem boundary (GIT_DISCOVERY_ACROSS_FILESYSTEM not set).
```
I tried different combination of SBATCH header as its seems it is related to the multi-threading and I get another very similar error message :
This is the error message for the job # is 37669225
```
[Error] [0:00:28] mpi-scheduler execution failed with error code 134
[Error] [0:00:28] Will now exit...
[Error] <class 'RuntimeError'> mpi-scheduler  execution failed with error code 134
Writing report file in /home/christec/Clematis/Run3/Analyses/Tree/raxml_all/report.txt
When reporting the issue, please always send us this file.
fatal: not a git repository (or any parent up to mount point /)
Stopping at filesystem boundary (GIT_DISCOVERY_ACROSS_FILESYSTEM not set).
```
I would be grateful if you could help me with this issue.
Best regards,
Camille Christe


## Post 2 by @Yann.Sagon (2025-01-28T15:55:08.923Z)

Dear @Camille.Christe[@Camille.Christe](https://hpc-community.unige.ch/u/camille.christe)
As I don’t know this software, I googled it and this is what I found: Parallelization · BenoitMorel/ParGenes Wiki · GitHub[Parallelization · BenoitMorel/ParGenes Wiki · GitHub](https://github.com/BenoitMorel/ParGenes/wiki/Parallelization) . Does it correspond to what you are talking about? If the
Can you please show your sbatch script?
Best
Yann


## Post 3 by @Yann.Sagon (2025-01-28T16:48:23.770Z)

Dear @Camille.Christe[@Camille.Christe](https://hpc-community.unige.ch/u/camille.christe)
As I don’t know this software, I googled it and this is what I found: Parallelization · BenoitMorel/ParGenes Wiki · GitHub[Parallelization · BenoitMorel/ParGenes Wiki · GitHub](https://github.com/BenoitMorel/ParGenes/wiki/Parallelization) . Does it correspond to what you are talking about?
Can you please show your sbatch script?
Best
Yann
