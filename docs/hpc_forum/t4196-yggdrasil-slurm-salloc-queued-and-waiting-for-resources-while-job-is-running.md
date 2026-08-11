# [yggdrasil][slurm] salloc - "queued and waiting for resources" while job is running

- Source: https://hpc-community.unige.ch/t/4196

- Created: 2026-01-16T14:53:09.559Z

- Tags: yggdrasil, slurm

- Posts: 1

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @maciej.falkiewicz (2026-01-16T14:53:09.610Z)

Dear @hpc[@hpc](https://hpc-community.unige.ch/u/hpc) Team,
My salloc job on Yggdrasil is still awaiting resources
```
(yggdrasil)-[falkiewi@login1 fixing-cfg]$ ./apptainer/salloc-cpu.sh $CPU_PARTITIONS 4:00:00 16000 4
salloc: Pending job allocation 43722373
salloc: job 43722373 queued and waiting for resources
```
while theoretically “running” at the same time:
```
(yggdrasil)-[falkiewi@login1 fixing-cfg]$ squeue -o "%.18i %.9P %.75j %.8u %.2t %.10M %.6D %.20R %q %S" | grep 43722373
          43722373 public-cp                                                                 interactive falkiewi  R       7:12      1               cpu027 normal 2026-01-16T15:44:59
```
This looks like a problem with slurm. Thank you in advance.
Kind regards,
Maciej
