# Partial cancellation of jobs, others completed fine

- Source: https://hpc-community.unige.ch/t/3496

- Created: 2024-06-17T19:57:14.977Z

- Posts: 1

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Mathieu.Dedenon (2024-06-17T19:57:15.018Z)

Hi HPC team,
## Primary informations
Username: dedenon
Cluster: Baobab
## Description
I ran a CPU job array last Friday on our private partition kruse-cpu, and I got cancellation of more than half of them, but the remaining ones are completed.
jobarray
jobarray851×1057 239 KB
[jobarray851×1057 239 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/b81e83470c2e3c1ad6bf78d23e5f3c9e31d07091.png)
Here is the sbatch file
sbatch
sbatch624×609 73.2 KB
[sbatch624×609 73.2 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/6102f0df183a2fb8b04bd2d5ca9cc4e2ddee798d.png)
I don’t think the issue is in my code, because I do parametric scanning for simulations with 20 independent realizations for 5 000 000 time steps, and some of them got completed (see slurm output below)
job 91
slurm91
slurm91923×398 93.4 KB
[slurm91923×398 93.4 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/d816db8cedaa1892661f020c4ef3d5eb7208774a.png)
but others are cancelled…
job 92
slurm92
slurm92949×333 84 KB
[slurm92949×333 84 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/d836bafc73530439c6a4320a54698eb0287c02a5.png)
Here are the seff reports for those jobs, this is neither TIMEOUT or OUT-OF-MEMORY issue
seff91-92
seff91-92541×530 78.5 KB
[seff91-92541×530 78.5 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/0f5a6746cbe2d74101fcb79c7d27d4690cd5096e.png)
Other example from job 115
slurm115
slurm115952×438 109 KB
[slurm115952×438 109 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/63b3d7d4595af7fbcce242559d0dddb7e422a3c9.png)
Finally I checked for node status in the partition, I got this
state-private-kruse
state-cpu195
## Steps to Reproduce
By essence it is not a kind of reproducible problem, this is something that already happened and I can’t figure out what triggers this… It happened several times and I just noticed that there is no issue if I launch small job arrays (~100), but large ones (~1000) get systematically cancelled.
One interesting information: looking at the slurm output files, all uncompleted jobs seem to have been cancelled at 08am on Saturday!
I will try on another partition to see if it changes anything.
Thanks by advance for your help on this !
Best,
Mathieu D.
