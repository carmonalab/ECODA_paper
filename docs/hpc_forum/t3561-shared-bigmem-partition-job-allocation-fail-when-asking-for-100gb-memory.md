# Shared-bigmem partition job allocation fail when asking for >100Gb memory

- Source: https://hpc-community.unige.ch/t/3561

- Created: 2024-07-23T08:40:48.860Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Felix.Hubert (2024-07-23T08:40:48.903Z)

Hello,
I am using the command
salloc --partition=“shared-bigmem” --mem=500000 --time=12:00:00
to ask for a job with a lot of memory (I need at least 400Gb), as I was used to do in the past year, but now it seems that there’s a threshold at 100Gb. The command
salloc --partition=“shared-bigmem” --mem=100000 --time=12:00:00
works, but for anything above 100Gb of memory, the waiting time for ressources is seemingly infinite.
Did you change the maximum threshold for the shared-bigmem partition? Or is it just temporary?
Thank you!
Félix


## Post 2 by @Adrien.Albert (2024-07-23T12:30:33.383Z)

Hi @Felix.Hubert[@Felix.Hubert](https://hpc-community.unige.ch/u/felix.hubert)
Felix.Hubert:
> salloc --partition=“shared-bigmem” --mem=500000 --time=12:00:00
- You have a syntax error the partition "shared-bigmem”
```
(baobab)-[alberta@login1 ~]$ salloc --partition=“shared-bigmem” --mem=500000 --time=12:00:00
salloc: error: invalid partition specified: “shared-bigmem”
salloc: error: Job submit/allocate failed: Invalid partition name specified
```
- It’s working for me:
```
(baobab)-[alberta@login1 ~]$ salloc --partition=shared-bigmem --mem=500000 --time=12:00:00
salloc: Pending job allocation 11506750
salloc: job 11506750 queued and waiting for resources
salloc: job 11506750 has been allocated resources
salloc: Granted job allocation 11506750
salloc: Nodes cpu186 are ready for job
(baobab)-[alberta@cpu186 ~]$
```
As you are required a large amount of memory you need to wait  until the resource is available and it’s depend on your priority:
https://doc.eresearch.unige.ch/hpc/slurm#job_priorities[https://doc.eresearch.unige.ch/hpc/slurm#job_priorities](https://doc.eresearch.unige.ch/hpc/slurm#job_priorities)
