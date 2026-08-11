# Cannot allocate cpu/gpu or submit any other job

- Source: https://hpc-community.unige.ch/t/3828

- Created: 2025-02-17T06:56:24.750Z

- Tags: baobab, slurm

- Posts: 5

- Category: 5

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Tomke.Schroeer (2025-02-17T06:56:24.806Z)

I tried to allocate a cpu and also tried multiple times to submit another job to baobab but every time it is cancelled right after submitting. The time limit does not exceed the time until the maintenance. I used the same scripts as before that used to work.
Did others experience the same and what could be the reason for it?
Thanks in advance
Tomke


## Post 2 by @Yann.Sagon (2025-02-17T07:39:48.001Z)

Dear @Tomke.Schroeer[@Tomke.Schroeer](https://hpc-community.unige.ch/u/tomke.schroeer)
please give us more details about the issue: on which cluster does it happens, an sbatch example, how you submit the sbatch if any.
hint: for the next time, if you use the category “support issues” on hpc-community you have a template to help you fill would help us help you.
Best


## Post 3 by @Tomke.Schroeer (2025-02-17T07:42:49.942Z)

Dear @Yann.Sagon[@Yann.Sagon](https://hpc-community.unige.ch/u/yann.sagon)
I am using baobab and for allocating a cpu I use this script:
```
#!/bin/bash
#SBATCH --job-name=bg-cpu-session
#SBATCH --partition=shared-cpu,private-dpnc-cpu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclude=cpu328
#SBATCH --time=10:00:00
#SBATCH --mem=64GB

echo "Allocated a background CPU session on node '`hostname -s`'."
srun sleep 12h
```


## Post 4 by @Yann.Sagon (2025-02-17T08:33:11.706Z)

Dear @Tomke.Schroeer[@Tomke.Schroeer](https://hpc-community.unige.ch/u/tomke.schroeer)
This is the reason:
```
(baobab)-[schroeer@login1 ~]$ beegfs-get-quota-home-scratch.sh
home dir: /home/users/s/schroeer
scratch dir: /srv/beegfs/scratch/users/s/schroeer

          user/group                 ||           size          ||    chunk files
  storage     |   name        |  id  ||    used    |    hard    ||  used   |  hard
  ----------------------------|------||------------|------------||---------|---------
home        |       schroeer|403388||    1.05 TiB| 1024.00 GiB||   663244|unlimited
scratch     |       schroeer|403388||   16.24 TiB|   unlimited||   624908| 10000000
(baobab)-[schroeer@login1 ~]$ touch aa
touch: cannot touch 'aa': Disk quota exceeded
```
By the way, to request  an interactive session on a compute node, it is better to use salloc. hpc:slurm [eResearch Doc][hpc:slurm [eResearch Doc]](https://doc.eresearch.unige.ch/hpc/slurm#interactive_jobs)


## Post 5 by @Tomke.Schroeer (2025-02-17T09:20:32.349Z)

Ok thanks, this solved it!
cheers
Tomke
