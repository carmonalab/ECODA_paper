# Gpu nodes not connecting to internet

- Source: https://hpc-community.unige.ch/t/3289

- Created: 2024-02-05T22:02:52.784Z

- Posts: 7

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nicola.Piasentin (2024-02-05T22:02:52.847Z)

Hello,
I was installing some modified version of GROMACS locally in my home folder. I salloc to a private-gervasio-gpu node and then module load the relevant modules and do the installation. I already did it successfully other times here on baobab so the problem is (relatively) new and wasn’t present before.
In brief, there is a passage where GROMACS has to downloads some libraries on its own. However, the passage fails because of time-out connection to the website where it retrieves the tarball. The website is up and working properly, so it’s not a website problem.
If I salloc to the node and try to wget something I get a timed-out connection, same for nodes in the shared-gpu partition. This doesn’t happen on the head node. As such, I think for some reasons the nodes are not able to connect via internet?
Best,
Nicola


## Post 2 by @Malte.Algren (2024-02-06T06:31:34.521Z)

I think I am finding the same issue on baobab. Some of us are using Weights&Biases to log jobs to the internet. I am on GPU044, and when I init wandb I get ConnectTimeout:
`wandb: Network error (ConnectTimeout), entering retry loop.`
`wandb: Network error (ConnectTimeout), entering retry loop.`
`Problem at: /home/users/a/algren/work/diffusion/run/run_diffusion.py 26 main `
`Error executing job with overrides: []`
`wandb: ERROR Run initialization has timed out after 90.0 sec.`


## Post 3 by @Adrien.Albert (2024-02-06T08:10:55.458Z)

Dear All,
Guilty!!! :pleading_face: I’ve made a modification on the admin node of each cluster that shouldn’t have any impact on external communication. It seems that this problem does not appear on Yggdrasil.
So I’ve undone the change to get everything working again:
```
(baobab)-[root@admin1 ~]$ clush -bw @compute fping google.ch

---------------
cpu[001-005,007-008,011,019,026-028,045-055,058-059,061-066,154-162,164-203,205-213,216-226,228-229,237-273,275-318,320-335],gpu[002,004-009,011-046] (242)
---------------
google.ch is alive
```
I apologize for the inconvenience :pray:t3:


## Post 4 by @Jingze.Duan (2024-02-06T08:59:14.766Z)

Hello,
I met the similar problem. My GROMACS jobs on yggdrasil all failed last night. And I resubmitted a job this morning. But I got the same error information as below.
Fatal error:
Cannot run short-ranged nonbonded interactions on a GPU because no GPU is
detected.
What should I do to fix it?
Best,
Jingze


## Post 5 by @Adrien.Albert (2024-02-06T10:52:28.827Z)

The fix has been applied this morning, you should not have this error anymore.


## Post 6 by @Jingze.Duan (2024-02-06T14:18:14.918Z)

Hi,
I resubmitted jobs this afternoon. One runs normally with gpu007, but another one only ran for a few sec on gpu008 and failed.
## Some lines of .out file:
Program:     gmx mdrun, version 2023.1
Source file: src/gromacs/taskassignment/findallgputasks.cpp (line 85)
MPI rank:    0 (out of 8)
Fatal error:
Cannot run short-ranged nonbonded interactions on a GPU because no GPU is
detected.
For more information and tips for troubleshooting, please check the GROMACS
website at Common Errors — GROMACS webpage https://www.gromacs.org documentation[Common Errors — GROMACS webpage https://www.gromacs.org documentation](http://www.gromacs.org/Documentation/Errors)
[1707228579.325814] [gpu008:1539109:0]           ib_md.c:1234 UCX  WARN  IB: ibv_fork_init() was disabled or failed, yet a fork() has been issued.
[1707228579.325821] [gpu008:1539109:0]           ib_md.c:1235 UCX  WARN  IB: data corruption might occur when using registered memory.
Best,
Jingze


## Post 7 by @Adrien.Albert (2024-02-06T14:36:51.702Z)

```
Fatal error:
Cannot run short-ranged nonbonded interactions on a GPU because no GPU is
detected.
```
This error does not seem to be link to network issue.
Please Could you , create a new thread and you share me your sbatch ?


## Post 8 by @Adrien.Albert (2024-02-06T14:37:49.699Z)
